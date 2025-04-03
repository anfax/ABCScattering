#!/bin/bash
# 构建并复制可执行文件到 ./run
fpm build --compiler ifort  --profile release && \
mkdir -p ./run && \
rm -f ./run/your_program && \
cp ./build/*/app/* ./run/