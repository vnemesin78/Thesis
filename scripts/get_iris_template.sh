#Config. files
PARAMETERS=`ls ../parametres_iris/*.cfg "../parametres_iris/image"/*.cfg`
APP_DIR=../Executables
$APP_DIR/get_iris_template.run ${1} ../Temp $PARAMETERS
