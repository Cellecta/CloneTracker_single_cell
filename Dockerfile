FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libglib2.0-0 \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md requirements.txt ./
COPY src ./src
COPY run_pipeline.py clonetracker_batch.py scRNAseq_QC.py ./
COPY best_sequence_umi.py barcode_process_umi_5.py barcode_final_assignment_50M.py ./

RUN pip install --upgrade pip && pip install .

ENTRYPOINT ["cellecta-full-pipeline"]
CMD ["--help"]
