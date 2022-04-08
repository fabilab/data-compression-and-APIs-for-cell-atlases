# Create service
aws lightsail create-container-service --service-name compressed-atlas --power small --scale 1

# Check service
aws lightsail get-container-services --service-name compressed-atlas

# Delete service
#aws lightsail delete-container-service --service-name compressed-atlas
