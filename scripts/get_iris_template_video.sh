#Config. files
PARAMETERS=`ls ../parametres_iris/*.cfg "../parametres_iris/video"/*.cfg`
APP_DIR=../Executables
$APP_DIR/get_iris_template_video.run ${1} ../Temp $PARAMETERS ${2}
