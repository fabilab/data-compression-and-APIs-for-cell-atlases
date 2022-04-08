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
# Flask app
COPY main.py .
# Backend modules
COPY helper.py .
# HTML files
COPY index.html .
COPY page2.html .
# Static assets
# FIXME: get rid of 10X data et similia
ADD static static

# Specify the command to run on container start
CMD [ "python", "./main.py" ]
