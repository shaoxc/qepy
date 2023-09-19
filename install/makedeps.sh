#!/bin/sh
LC_ALL=C
export LC_ALL
QEDIR=${qedir}
if [ ! $QEDIR ];then
	QEDIR=../../../
fi
#QESUBS=`ls -d ${QEDIR}/*/`
${QEDIR}/install/moduledep.sh  > make.depend
${QEDIR}/install/includedep.sh >> make.depend

# handle special cases
sed -i '/@/d'                 make.depend
grep -v -x -f make.depend.cetddft make.depend > make.depend.tmp
mv make.depend.tmp make.depend
