# UbE3 Ligase Activity Profiling Analysis
The Python package can be installed through executing the following command in a Python console:

    python3 -m pip install ube3_apa

After successful installation, you may test the code with the [testing data](https://github.com/Chenlab-UMN/Ub-E3-ligase-Activity-Profiling-Analysis/tree/master/testdata)

# Getting Started
Each script will begin with
    
    import ube3_apa
    input_directory = '<directory_of_testdata_folder>'

Next, there are different ways the analysis can be run.

## Standard Analysis

    ube3_apa.e3enrich(siteratio_dir = input_directory + '/testdata/siteratio_testdata.csv', input_type = 'UniprotAC', output_dir = '<desired_output_directory>', exp_label = '1', grouped = False, output_ratio = True, log2trans = False)

## Grouped Analysis

    ube3_apa.e3enrich(siteratio_dir = input_directory + '/testdata/siteratio_testdata.csv', input_type = 'UniprotAC', output_dir = '<desired_output_directory>', exp_label = '2', grouped = True, output_ratio = True, log2trans = False)

## Analysis with Protein Normalization

    ube3_apa.e3enrich(siteratio_dir = input_directory + '/testdata/siteratio_testdata.csv', input_type = 'UniprotAC', output_dir = '<desired_output_directory>', exp_label = '3', proratio_dir = input_directory + '/testdata/proteinratio_testdata.csv', grouped = False, output_ratio = True, log2trans = False)

After running the code above, you will find some csv files in the output directory. The enrichment p-values are listed in files initiated with “UbE3_APA”, and there are also other related data such as the number of substrates found and average site ratios included.

# Parameters
## General Command
    e3enrich(siteratio_dir, input_type, output_dir, exp_label='', proratio_dir="None", grouped=False, output_ratio=False, log2trans=True)
Perform Ub E3 ligase activity profiling analysis based on ratio of E3 ligase substrates Parameters can be customized through key-value pairs as below.
| Parameter | Type | Description |
| --- | --- | --- |
| `siteratio_dir` | String | The directory of the file that contains the ratio of every site. |
| `input_type` | String \in {'UniprotAC', 'protein', 'gene symbol', 'gene'} | The type of IDs used in the site ratio file and the protein ratio file, format in both files should be the same. <br><br>"UniprotAC" or "protein" example: Q00987, P40337, Q9HAU4, O43791 <br><br>"gene symbol" or "gene" example: MDM2, VHL, SMUF2, SPOP List of valid input and examples will be shown if it is an invalid value. |
| `output_dir` | String | The directory where E3 enrichment result files will be generated. |
| `exp_label` | String, default: '' | The string attached to the output file name that separates different results when there are multiple groups. |
| `proratio_dir` | String, default: '' | The directory of the file that contains the ratio of every protein. A valid directory input here will trigger normalization by corresponding protein ratio for all files in this run. By default, the output will not be normalized by protein ratio. |
| `grouped` | Boolean, default: False | If True, enrichment results will show leading E3 ligase and grouped E3 ligase in each row of the result instead of showing each E3 ligase individually. E3 ligases are grouped according to the relationship of the detected substrate in this run. |
| `ratio_output` | Boolean, default: False | If True, files that contain the ratio of each ubiquitylation site and the average ratio of each ubiquitinated protein will be generated. |
| `log2trans` | Boolean, default: True | If True, ratios will take transform y=log2(x) before the E3 ligase enrichment analysis. If False, the ratio from files will be used for enrichment analysis directly. |
# Output
This function generates files based on E3 ligase enrichment results and does not have any returned values.
# Reference
**UbE3-APA: A Bioinformatic Strategy to Elucidate Ubiquitin E3 Ligase Activities in Quantitative Proteomics Study**
<br>Yao Gong, Yue Chen
<br>bioRxiv 2022.01.16.476541; doi: https://doi.org/10.1101/2022.01.16.476541