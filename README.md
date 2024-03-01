# scRNAseq
**Author**: Daelin Jensen  
**Date**: Spring 2024

#### Description:

[Seurat Documentation](https://satijalab.org/seurat/)

#### Set-up:
1. Clone git repository locally  
`git clone https://github.com/Jensen416/scRNAseq.git`
2. Create metadata file  
It is imperative that the metadata file is formatted correctly. See below:  
| file_path               | sample_name | age  | sex  | treatment |  
|---|---|---|---|---|  
| C:/usr/dir/sample_01.h5 | sample_01 | aged | male | control |  
| C:/usr/dir/sample_02.h5 | sample_02 | aged | male | drug |
| ... | ... | ... | ... | ... |  

| Header 1 | Header 2 | Header 3 |
|---|---|---|
| Content 1, row 1 | Content 2, row 1 | Content 3, row 1 |
| Content 1, row 2 | Content 2, row 2 | Content 3, row 2 |
| ... | ... | ... |
