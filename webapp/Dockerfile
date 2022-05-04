# Set base image (host OS)
FROM python:3.10

# By default, listen on port 5000
EXPOSE 5000/tcp

# Set the working directory in the container
WORKDIR /app

# Copy the dependencies file to the working directory
COPY requirements.txt .

# Install any dependencies
RUN pip install -r requirements.txt

# Copy the content of the local src directory to the working directory
# Python/Flask backend
COPY app.py .
COPY ca_API.py .
ADD models models
ADD text_recognition text_recognition
ADD validation validation
ADD voice_recognition voice_recognition

# Templates
ADD templates templates

# Static assets
# FIXME: get rid of 10X data et similia
ADD static static

# Specify the command to run on container start
CMD [ "python", "./app.py" ]
