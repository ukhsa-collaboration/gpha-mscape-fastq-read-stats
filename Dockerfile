FROM python:3.12

ADD . /app
WORKDIR /app
RUN pip install .

CMD ["fastq_read_stats"]
