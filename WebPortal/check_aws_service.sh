# Create service
aws lightsail create-container-service --service-name compressed-atlas-service --power small --scale 1

# Check service
aws lightsail get-container-services --service-name compressed-atlas-service

# Push image to service
aws lightsail push-container-image --service-name compressed-atlas-service --label compressed-atlas --image compressed-atlas

# Deploy container onto a public URL
aws lightsail create-container-service-deployment --service-name compressed-atlas-service --containers file://containers.json --public-endpoint file://public-endpoint.json

# Delete service
#aws lightsail delete-container-service --service-name compressed-atlas
