Bootstrap: docker
from: condaforge/mambaforge

%labels
    MAINTAINER jallen17@illinois.edu

%files
    /home/joshfactorial/NEAT

%post
    cd /home/joshfactorial/NEAT/ && \
    apt update -y && apt upgrade -y &&\
    apt install gcc -y && apt install g++ -y && \
    mamba env create -f environment.yml && \
    . /opt/conda/etc/profile.d/conda.sh && \
    conda activate neat && \
    poetry build -f wheel && \
    pip install dist/*whl && \
    rm -rf neat/dist && \
    rm -rf /root/.cache/pypoetry && \
    pip cache purge && \
    conda clean -afy && \
    mkdir /data

%environment
    export PATH=/opt/conda/envs/neat/bin:$PATH

%runscript
    echo "Running neat"
    neat gen_frag_model \
        -i /home/jallen17/neat_data/hg38_aligned.sorted.dedupped.bam \
        -o /home/jallen17/testing_frags/new_fraglen \
        --min_reads 100
    echo "neat finished"