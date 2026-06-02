ARG ROOT_IMAGE=rootproject/root:6.38.00-ubuntu25.10
FROM ${ROOT_IMAGE}

WORKDIR /workspace

COPY requirements.txt /tmp/requirements.txt

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       software-properties-common \
       python3-pip \
       python3-dev \
       python3-venv \
    && add-apt-repository universe \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
       xrootd-client \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv --system-site-packages /opt/venv \
    && /opt/venv/bin/python -m pip install --no-cache-dir --upgrade pip \
    && /opt/venv/bin/python -m pip install --no-cache-dir -r /tmp/requirements.txt

COPY . /workspace

ENV PATH="/opt/venv/bin:${PATH}"
ENV PYTHONPATH=/workspace/src

COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

CMD ["bash"]
