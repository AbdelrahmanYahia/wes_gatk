mkdir -p empty_sample_names
cat sample_sheet.txt | while read line ; do  touch empty_sample_names/$line; done