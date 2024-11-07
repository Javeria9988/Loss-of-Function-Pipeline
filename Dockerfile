# Start from the official BAKTA image
FROM oschwengers/bakta:latest

# Set environment variable for the database location
ENV BAKTA_DB=/db-light

# Copy the local db-light directory into the image
COPY db-light /db-light

# Optional: Set a default command (you can modify this based on usage)
CMD ["bash"]
