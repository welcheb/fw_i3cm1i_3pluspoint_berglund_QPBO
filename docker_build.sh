# build docker image using Dockerfile in current folder
USER="welcheb"
REPO="fw_i3cm1i_3pluspoint_berglund_qpbo"
NAME="$USER/$REPO"
docker build -t=$NAME .

# save tar.gz
#docker save $NAME | gzip -9 > ./docker_save/$REPO.tar.gz
