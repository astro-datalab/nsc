#!/bin/bash
if [ "$#" -eq 1 ]; then
  export FILE=${1}
else
  export FILE="archive_data"
fi
echo "Writing NOAO Archive search results to >>${FILE}<<"
\rm ${FILE}.txt ${FILE}.fits ${FILE}.fits.gz >& /dev/null
echo "Input password for dbreader:"
read -s password
export PGPASSWORD=${password}
export PSQL="psql -h db.sdm.noao.edu -p 5432 -U dbreader -d metadata"
# Get header line
$PSQL -E  -c "select b.instrument,a.sb_id,a.sb_recno,a.dtnsanam,a.dtacqnam,a.data_product_id, \
                a.uri,b.prop_id,b.ra,b.dec,b.exposure,b.release_date,b.date_obs,b.filter, \
                b.mjd_obs,b.obstype,b.plver,b.proctype,b.prodtype,b.filename,b.pldname \
                from edu_noao_nsa.data_product a, voi.siap b where a.data_product_id=b.fits_data_product_id limit 1;" | head -1  > ${FILE}.txt
# Get the data
$PSQL -E -t -c "select b.instrument,a.sb_id,a.sb_recno,a.dtnsanam,a.dtacqnam,a.data_product_id, \
                a.uri,b.prop_id,b.ra,b.dec,b.exposure,b.release_date,b.date_obs,b.filter, \
                b.mjd_obs,b.obstype,b.plver,b.proctype,b.prodtype,b.filename,b.pldname \
                from edu_noao_nsa.data_product a, voi.siap b where a.data_product_id=b.fits_data_product_id;" >> ${FILE}.txt
if [ -s ${FILE}.txt ]; then
  python -c "from astropy.table import Table; dat=Table.read('${FILE}.txt',format='ascii',delimiter='|'); dat.write('${FILE}.fits')"
  gzip ${FILE}.fits
  \rm ${FILE}.txt
else
  echo "No ${FILE}.txt file"
fi
