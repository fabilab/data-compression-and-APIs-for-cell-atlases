# Set base image (host OS)
FROM python:3.10

# By default, listen on port 5000
EXPOSE 5000/tcp

# Set the working directory in the container
WORKDIR /app

# Copy the dependencies file to the working directory
COPY requirements.txt .
COPY secret_key.txt .

# Install any dependencies
RUN pip install -r requirements.txt

# Templates
ADD templates templates

# Static assets
ADD static static

# Modules for backend
ADD api api
ADD models models
ADD text_recognition text_recognition
ADD validation validation
ADD voice_recognition voice_recognition

# Google cloud key (NOTE: careful!!)
# NOTE: CAREFUL CAREFUL
ADD gcloud_key gcloud_key
# NOTE: CAREFUL CAREFUL

# App file
COPY app.py .

# Specify the command to run on container start
CMD [ "python", "./app.py" ]
