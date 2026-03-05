import json
import logging
from rabbit_mq.rabbit_mq_consumer import create_consumer
from services.csv_top_10_rows_service import get_top_10_rows_from_output
from services.download_service import download_file
from services.message_producer import create_producer
from pipelines import run_pipeline
from services.s3_upload_service import process_and_upload_analysis
from services.clear_service import clear_all_folders

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load configuration
with open("config.json", "r") as f:
    config = json.load(f)


def process_request(message_data):
    """
    Process the consumed message.

    Args:
        message_data (dict): The jsonified message data
    """
    logger.info(f"Processing request: {message_data}")
    logger.info(f"Invoking the Download Service for {message_data['AnalysisID']}")

    # Create message producer instance
    producer = create_producer(config)

    try:
        if "expression_file_url" in message_data and "mapping_file_url" in message_data:
            # Download Files
            expression_path = download_file(message_data["expression_file_url"])
            mapping_path = download_file(message_data["mapping_file_url"])
            # Invoke Pipeline
            run_pipeline.main(
                expression_file=expression_path,
                mapping_file=mapping_path,
                config_path="./config.json",
            )

            # Upload Files to Cloudflare R2 and retrieve urls
            response_data = process_and_upload_analysis(message_data["AnalysisID"])

            # response_data = get_top_10_rows_from_output("Output", dest=response_data)
            # Print the formatted JSON output
            print("\nFinal Response Message:")
            print(json.dumps(response_data, indent=4))
            # Publish response message to response_queue
            correlation_id = message_data.get("AnalysisID", None)
            success = producer.publish_message(
                response_data, correlation_id=correlation_id
            )

            if success:
                logger.info(
                    f"Response published successfully for Analysis ID: {correlation_id}"
                )
            else:
                logger.error(
                    f"Failed to publish response for Analysis ID: {correlation_id}"
                )
            clear_all_folders()
        else:
            error_response = {
                "AnalysisID": message_data.get("AnalysisID", "unknown"),
                "status": "error",
                "message": "Missing required files: expression_file_url or mapping_file_url",
            }
            producer.publish_message(
                error_response, correlation_id=message_data.get("AnalysisID")
            )
            logger.error(f"Missing required file URLs in request: {message_data}")

    except Exception as e:
        # Send error response to response_queue
        error_response = {
            "AnalysisID": message_data.get("AnalysisID", "unknown"),
            "status": "error",
            "message": str(e),
        }
        producer.publish_message(
            error_response, correlation_id=message_data.get("AnalysisID")
        )
        logger.error(f"Error processing request: {e}")

    finally:
        producer.close()


if __name__ == "__main__":
    # Create consumer instance
    consumer = create_consumer(config)

    # Start consuming messages from request_queue and jsonify them
    try:
        consumer.consume(callback=process_request)
    except Exception as e:
        logger.error(f"Consumer error: {e}")
    finally:
        consumer.close()
