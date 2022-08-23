Some of the files were zipped due to github file size restriction.   
The file rnd1_features.csv was split into 3, and can be retrieved using the python code below.

``` python
import pandas as pd  

one = pd.read_csv('rnd1_features__PART1.csv.gz')
two = pd.read_csv('rnd1_features__PART2.csv.gz') 
three = pd.read_csv('rnd1_features__PART3.csv.gz')
features = pd.concat([one, two, three])
features.to_csv('rnd1_features.csv')
```
