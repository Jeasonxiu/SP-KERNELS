#!/bin/bash
#
#ln -sf ../../../../KERNEL.BAK_9_30_2016/OUTPUT_FILES_200000-0.01-2500-23 OUTPUTS_BAK
rm OUTPUTS_BAK OUTPUTS_CUR
ln -sf ../OUTPUT_FILES_100000-0.01-7500-15-DBLPERIOD-SALVAGE OUTPUTS_BAK
ln -sf ../OUTPUT_FILES_100000-0.01-7500-15-DBLPERIOD OUTPUTS_CUR
#
echo
echo Par_file.in
echo
diff OUTPUTS_BAK/Par_file.in OUTPUTS_CUR/Par_file.in
echo
echo
echo SOURCE
echo
diff OUTPUTS_BAK/DATA/SOURCE OUTPUTS_CUR/DATA/SOURCE
echo
echo
echo interfaces
echo
diff OUTPUTS_BAK/interfaces.dat OUTPUTS_CUR/interfaces.dat
echo
echo
echo tomofile
echo
diff OUTPUTS_BAK/DATA/tomo_file.xyz OUTPUTS_CUR/DATA/tomo_file.xyz
