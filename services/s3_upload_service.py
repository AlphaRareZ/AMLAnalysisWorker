import os
import json
import boto3
from botocore.config import Config

# Load config
with open('config.json', 'r') as f:
    config = json.load(f)

ACCOUNT_ID = config["CFR2"]["ACCOUNT_ID"]
ACCESS_KEY = config["CFR2"]["ACCESS_KEY"]
SECRET_KEY = config["CFR2"]["SECRET_KEY"]
BUCKET_NAME = config["CFR2"]["BUCKET_NAME"]
PUBLIC_DOMAIN = config["CFR2"]["PUBLIC_DOMAIN"] if "PUBLIC_DOMAIN" in config["CFR2"] else None
# Initialize S3 Client
s3 = boto3.client(
    service_name="s3",
    endpoint_url=f"https://{ACCOUNT_ID}.r2.cloudflarestorage.com",
    aws_access_key_id=ACCESS_KEY,
    aws_secret_access_key=SECRET_KEY,
    region_name="auto",
    config=Config(signature_version="s3v4"),
)

def process_and_upload_analysis(analysis_id):
    """
    Navigates through the local 'Output' folder, uploads each file to a 
    folder named after the Analysis ID in the bucket, and returns their URLs.
    """
    local_folder = "Output" # You can change this if your local path includes the ID, e.g., f"Output/{analysis_id}"
    
    if not os.path.exists(local_folder):
        print(f"Error: Local folder '{local_folder}' does not exist.")
        return {"success": False, "analysis_id": analysis_id, "links": {}, "error": "Output folder not found"}

    result_links = {}
    uploaded_count = 0
    failed_count = 0

    # 1. Navigate through the Output folder
    for root, _, files in os.walk(local_folder):
        for file in files:
            local_path = os.path.join(root, file)
            
            # Preserve folder structure relative to the 'Output' folder
            relative_path = os.path.relpath(local_path, local_folder)
            relative_path = relative_path.replace("\\", "/") # Windows compatibility
            
            # 2. Make a folder in the bucket for the analysis ID
            object_name = f"uploads/{analysis_id}/{relative_path}"
            
            try:
                # Upload the file
                s3.upload_file(local_path, BUCKET_NAME, object_name)
                print(f"✓ Uploaded: {object_name}")
                uploaded_count += 1
                
                # Check file extension to determine expiration strategy
                file_ext = os.path.splitext(file)[1].lower()
                
                if file_ext == '.csv':
                    # CSV: 1 Hour Expiration
                    url = s3.generate_presigned_url(
                        ClientMethod="get_object",
                        Params={"Bucket": BUCKET_NAME, "Key": object_name},
                        ExpiresIn=3600, 
                    )
                elif file_ext in ['.png', '.jpg', '.jpeg', '.gif']:
                    # Images: "Unlimited"
                    # OPTION A: Max allowed limit for presigned URLs (7 days)
                    url = s3.generate_presigned_url(
                        ClientMethod="get_object",
                        Params={"Bucket": BUCKET_NAME, "Key": object_name},
                        ExpiresIn=604800, # 7 days (604,800 seconds)
                    )
                    
                    # OPTION B: Truly unlimited (Uncomment if you have a Public R2 Custom Domain set up)
                    # PUBLIC_DOMAIN = "https://pub-yourdomain.r2.dev"
                    # url = f"{PUBLIC_DOMAIN}/{object_name}"
                else:
                    # Default expiration for other types (e.g., 24 hours)
                    url = s3.generate_presigned_url(
                        ClientMethod="get_object",
                        Params={"Bucket": BUCKET_NAME, "Key": object_name},
                        ExpiresIn=86400,
                    )
                
                # Add the URL to our response dictionary using the filename as the key
                result_links[relative_path] = url
                
            except Exception as e:
                print(f"✗ Failed to upload {local_path}: {str(e)}")
                failed_count += 1
                result_links[relative_path] = f"Error: {str(e)}"

    # Build the final response dictionary
    return {
        "success": failed_count == 0 and uploaded_count > 0,
        "analysis_id": analysis_id,
        "files_uploaded": uploaded_count,
        "files_failed": failed_count,
        "download_links": result_links
    }

# ==========================================
# Example Usage:
# ==========================================
if __name__ == "__main__":
    # Just pass the Analysis ID!
    response_data = process_and_upload_analysis("REQ-2024-001")
    
    # Print the formatted JSON output (which you can pass straight to RabbitMQ)
    print("\nFinal Response Message:")
    print(json.dumps(response_data, indent=4))