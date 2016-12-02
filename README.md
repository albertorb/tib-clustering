

# Smart Techniques in Bioinformatics

#### In this notebook, we explore the data and face MA against clustering to identify those proteins whose presence or absence could affect on human infarcts.

First, we use MA technique to capture outliers. That give us an aproximation of what we are searching for. Then we use k-means clustering technique trying to achieve more accurate results.

### Data loading


```R
dts <- read.csv('hg19.cage_peak_phase1and2combined_tpm_ann_decoded.osc.txt.gz.extract.tsv',sep = '\t' ,header=TRUE)
head(dts)
```


<table>
<thead><tr><th scope=col>X00Annotation</th><th scope=col>short_description</th><th scope=col>uniprot_id</th><th scope=col>Astrocyte...cerebellum..donor1.CNhs11321.11500.119F6</th><th scope=col>Astrocyte...cerebral.cortex..donor1.CNhs10864.11235.116D2</th><th scope=col>brain..adult..donor1.CNhs11796.10084.102B3</th><th scope=col>brain..adult..pool1.CNhs10617.10012.101C3</th><th scope=col>brain..fetal..pool1.CNhs11797.10085.102B4</th><th scope=col>breast..adult..donor1.CNhs11792.10080.102A8</th><th scope=col>cerebellum...adult..donor10196.CNhs13799.10173.103C2</th><th scope=col>⋯</th><th scope=col>thalamus...adult..donor10196.CNhs13794.10168.103B6</th><th scope=col>thalamus..adult..donor10252.CNhs12314.10154.103A1</th><th scope=col>thalamus..adult..donor10258..tech_rep1.CNhs14223.10370.105G1</th><th scope=col>thalamus..adult..donor10258..tech_rep2.CNhs14551.10370.105G1</th><th scope=col>thalamus..newborn..donor10223.CNhs14084.10366.105F6</th><th scope=col>throat..fetal..donor1.CNhs11770.10061.101H7</th><th scope=col>thyroid..fetal..donor1.CNhs11769.10060.101H6</th><th scope=col>tongue.epidermis..fungiform.papillae...donor1.CNhs13460.10288.104F9</th><th scope=col>umbilical.cord..fetal..donor1.CNhs11765.10057.101H3</th><th scope=col>uterus..fetal..donor1.CNhs11763.10055.101H1</th></tr></thead>
<tbody>
	<tr><td>chr10:100013403..100013414,-  </td><td>p@chr10:100013403..100013414,-</td><td>NA                            </td><td>0.00                          </td><td> 0.00                         </td><td>0                             </td><td>0                             </td><td>0                             </td><td>0.00                          </td><td>0                             </td><td>⋯                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td></tr>
	<tr><td>chr10:100027943..100027958,-  </td><td>p1@LOXL4                      </td><td>uniprot:Q96JB6                </td><td>0.12                          </td><td>11.45                         </td><td>0                             </td><td>0                             </td><td>0                             </td><td>2.17                          </td><td>0                             </td><td>⋯                             </td><td>5.80                          </td><td>0.31                          </td><td>5.65                          </td><td>2.99                          </td><td>0                             </td><td>1.19                          </td><td>1.01                          </td><td>2.58                          </td><td>7.04                          </td><td>4.48                          </td></tr>
	<tr><td>chr10:100076685..100076699,+  </td><td>p@chr10:100076685..100076699,+</td><td>NA                            </td><td>0.00                          </td><td> 0.00                         </td><td>0                             </td><td>0                             </td><td>0                             </td><td>0.00                          </td><td>0                             </td><td>⋯                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td></tr>
	<tr><td>chr10:100150910..100150935,-  </td><td>p@chr10:100150910..100150935,-</td><td>NA                            </td><td>0.00                          </td><td> 0.00                         </td><td>0                             </td><td>0                             </td><td>0                             </td><td>0.00                          </td><td>0                             </td><td>⋯                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td></tr>
	<tr><td>chr10:100150951..100150962,-  </td><td>p@chr10:100150951..100150962,-</td><td>NA                            </td><td>0.00                          </td><td> 0.00                         </td><td>0                             </td><td>0                             </td><td>0                             </td><td>0.00                          </td><td>0                             </td><td>⋯                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td></tr>
	<tr><td>chr10:100150986..100150988,+  </td><td>p@chr10:100150986..100150988,+</td><td>NA                            </td><td>0.00                          </td><td> 0.00                         </td><td>0                             </td><td>0                             </td><td>0                             </td><td>0.00                          </td><td>0                             </td><td>⋯                             </td><td>0.83                          </td><td>0.00                          </td><td>0.00                          </td><td>0.64                          </td><td>0                             </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td><td>0.00                          </td></tr>
</tbody>
</table>



### Data cleaning

#### Shorten info

We remove the first and the second row because it is enough just with uniprot_id for our purpose.


```R
dts <- dts[,-1:-2]
cat("Number of gens: ", nrow(dts), "\n") 
cat("Number of tissues: ", ncol(dts)) 
```

    Number of gens:  59173 
    Number of tissues:  71

#### Remove missing values


```R
dts <- na.omit(dts)
dts <- subset(dts, uniprot_id!="")
head(dts)
```


<table>
<thead><tr><th></th><th scope=col>uniprot_id</th><th scope=col>Astrocyte...cerebellum..donor1.CNhs11321.11500.119F6</th><th scope=col>Astrocyte...cerebral.cortex..donor1.CNhs10864.11235.116D2</th><th scope=col>brain..adult..donor1.CNhs11796.10084.102B3</th><th scope=col>brain..adult..pool1.CNhs10617.10012.101C3</th><th scope=col>brain..fetal..pool1.CNhs11797.10085.102B4</th><th scope=col>breast..adult..donor1.CNhs11792.10080.102A8</th><th scope=col>cerebellum...adult..donor10196.CNhs13799.10173.103C2</th><th scope=col>cerebellum..adult..donor10252.CNhs12323.10166.103B4</th><th scope=col>cerebellum..newborn..donor10223.CNhs14075.10357.105E6</th><th scope=col>⋯</th><th scope=col>thalamus...adult..donor10196.CNhs13794.10168.103B6</th><th scope=col>thalamus..adult..donor10252.CNhs12314.10154.103A1</th><th scope=col>thalamus..adult..donor10258..tech_rep1.CNhs14223.10370.105G1</th><th scope=col>thalamus..adult..donor10258..tech_rep2.CNhs14551.10370.105G1</th><th scope=col>thalamus..newborn..donor10223.CNhs14084.10366.105F6</th><th scope=col>throat..fetal..donor1.CNhs11770.10061.101H7</th><th scope=col>thyroid..fetal..donor1.CNhs11769.10060.101H6</th><th scope=col>tongue.epidermis..fungiform.papillae...donor1.CNhs13460.10288.104F9</th><th scope=col>umbilical.cord..fetal..donor1.CNhs11765.10057.101H3</th><th scope=col>uterus..fetal..donor1.CNhs11763.10055.101H1</th></tr></thead>
<tbody>
	<tr><th scope=row>2</th><td>uniprot:Q96JB6                              </td><td> 0.12                                       </td><td>11.45                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 2.17                                       </td><td> 0.00                                       </td><td> 0.22                                       </td><td> 1.03                                       </td><td>⋯                                           </td><td> 5.80                                       </td><td> 0.31                                       </td><td> 5.65                                       </td><td> 2.99                                       </td><td> 0.00                                       </td><td> 1.19                                       </td><td> 1.01                                       </td><td> 2.58                                       </td><td> 7.04                                       </td><td> 4.48                                       </td></tr>
	<tr><th scope=row>7</th><td>uniprot:Q8N2H3                              </td><td> 7.51                                       </td><td> 6.30                                       </td><td> 3.88                                       </td><td> 3.71                                       </td><td> 2.00                                       </td><td> 5.07                                       </td><td> 1.53                                       </td><td> 1.99                                       </td><td> 6.72                                       </td><td>⋯                                           </td><td> 4.15                                       </td><td> 7.34                                       </td><td> 4.23                                       </td><td> 3.84                                       </td><td> 4.28                                       </td><td>15.24                                       </td><td>11.92                                       </td><td> 7.74                                       </td><td> 7.04                                       </td><td>15.87                                       </td></tr>
	<tr><th scope=row>8</th><td>uniprot:Q8N2H3                              </td><td> 3.69                                       </td><td> 3.43                                       </td><td> 1.94                                       </td><td> 0.65                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 1.53                                       </td><td> 0.55                                       </td><td> 1.55                                       </td><td>⋯                                           </td><td> 0.83                                       </td><td> 2.69                                       </td><td> 2.82                                       </td><td> 1.92                                       </td><td> 0.00                                       </td><td> 3.33                                       </td><td> 3.03                                       </td><td> 2.58                                       </td><td> 2.35                                       </td><td> 5.29                                       </td></tr>
	<tr><th scope=row>15</th><td>uniprot:Q92902,uniprot:Q658M9,uniprot:Q8WXE5</td><td>35.58                                       </td><td>20.61                                       </td><td>30.05                                       </td><td>21.71                                       </td><td>19.96                                       </td><td>31.90                                       </td><td>13.73                                       </td><td>20.79                                       </td><td>21.72                                       </td><td>⋯                                           </td><td>15.75                                       </td><td>26.89                                       </td><td>24.00                                       </td><td>26.24                                       </td><td>12.83                                       </td><td>19.76                                       </td><td>22.03                                       </td><td>12.90                                       </td><td>32.85                                       </td><td>17.09                                       </td></tr>
	<tr><th scope=row>24</th><td>uniprot:Q8WWQ2                              </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 1.02                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 2.99                                       </td><td> 0.00                                       </td><td>⋯                                           </td><td> 5.80                                       </td><td> 8.38                                       </td><td> 4.23                                       </td><td> 8.96                                       </td><td> 0.71                                       </td><td> 0.71                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 7.04                                       </td><td> 5.70                                       </td></tr>
	<tr><th scope=row>25</th><td>uniprot:Q8WWQ2                              </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td>⋯                                           </td><td> 1.66                                       </td><td> 0.10                                       </td><td> 0.00                                       </td><td> 0.43                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td><td> 0.00                                       </td></tr>
</tbody>
</table>



##### Group repeated experiments of the same protein

As we can see there are some repeated experiments so we group them and store the mean of their gene expresion.
There are some experiments in which more than one gen was measured but we will treat them as if they were a single gen.


```R
library(plyr)
small_dts <- ddply(dts, 'uniprot_id', numcolwise(mean))

```

#### Change index

Now, to work more confortable, we assign the protein's id as indexes of the data frame. In addition, that ensure that experiments are unique and therefore step above worked.


```R
rownames(small_dts) <- small_dts$uniprot_id
```

uniprot_id is then removed from columns


```R
small_dts <- small_dts[,-1]
tail(small_dts)
```


<table>
<thead><tr><th></th><th scope=col>Astrocyte...cerebellum..donor1.CNhs11321.11500.119F6</th><th scope=col>Astrocyte...cerebral.cortex..donor1.CNhs10864.11235.116D2</th><th scope=col>brain..adult..donor1.CNhs11796.10084.102B3</th><th scope=col>brain..adult..pool1.CNhs10617.10012.101C3</th><th scope=col>brain..fetal..pool1.CNhs11797.10085.102B4</th><th scope=col>breast..adult..donor1.CNhs11792.10080.102A8</th><th scope=col>cerebellum...adult..donor10196.CNhs13799.10173.103C2</th><th scope=col>cerebellum..adult..donor10252.CNhs12323.10166.103B4</th><th scope=col>cerebellum..newborn..donor10223.CNhs14075.10357.105E6</th><th scope=col>dura.mater..adult..donor1.CNhs10648.10041.101F5</th><th scope=col>⋯</th><th scope=col>thalamus...adult..donor10196.CNhs13794.10168.103B6</th><th scope=col>thalamus..adult..donor10252.CNhs12314.10154.103A1</th><th scope=col>thalamus..adult..donor10258..tech_rep1.CNhs14223.10370.105G1</th><th scope=col>thalamus..adult..donor10258..tech_rep2.CNhs14551.10370.105G1</th><th scope=col>thalamus..newborn..donor10223.CNhs14084.10366.105F6</th><th scope=col>throat..fetal..donor1.CNhs11770.10061.101H7</th><th scope=col>thyroid..fetal..donor1.CNhs11769.10060.101H6</th><th scope=col>tongue.epidermis..fungiform.papillae...donor1.CNhs13460.10288.104F9</th><th scope=col>umbilical.cord..fetal..donor1.CNhs11765.10057.101H3</th><th scope=col>uterus..fetal..donor1.CNhs11763.10055.101H1</th></tr></thead>
<tbody>
	<tr><th scope=row>uniprot:Q9Y6Y8,uniprot:F5H0L8</th><td>29.55 </td><td>34.63 </td><td>35.87 </td><td>26.530</td><td>31.94 </td><td>42.77 </td><td>64.08 </td><td>34.17 </td><td>42.93 </td><td>43.730</td><td>⋯     </td><td>53.07 </td><td>24.72 </td><td>31.06 </td><td>41.81 </td><td>43.47 </td><td>51.66 </td><td>62.450</td><td>51.60 </td><td>44.58 </td><td>39.060</td></tr>
	<tr><th scope=row>uniprot:Q9Y6Y9,uniprot:D0VAW8</th><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.185</td><td> 0.00 </td><td> 0.72 </td><td> 0.00 </td><td> 0.11 </td><td> 0.26 </td><td> 1.485</td><td>⋯     </td><td> 0.83 </td><td> 0.52 </td><td> 0.00 </td><td> 0.00 </td><td> 1.07 </td><td> 0.12 </td><td> 0.305</td><td> 0.00 </td><td> 0.00 </td><td> 0.405</td></tr>
	<tr><th scope=row>uniprot:Q9Y6Y9,uniprot:E5RJJ7</th><td> 0.49 </td><td> 0.29 </td><td> 2.91 </td><td>13.080</td><td> 2.00 </td><td>18.12 </td><td> 1.53 </td><td> 3.10 </td><td> 4.66 </td><td>29.190</td><td>⋯     </td><td>11.61 </td><td>12.10 </td><td> 7.06 </td><td> 7.25 </td><td>14.96 </td><td> 5.24 </td><td> 1.410</td><td> 2.58 </td><td> 9.39 </td><td> 7.730</td></tr>
	<tr><th scope=row>uniprot:Q9Y6Z2</th><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td><td>⋯     </td><td> 0.00 </td><td> 0.21 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 1.67 </td><td> 3.230</td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td></tr>
	<tr><th scope=row>uniprot:Q9Y6Z5</th><td> 0.37 </td><td> 0.86 </td><td> 0.00 </td><td> 0.280</td><td> 9.98 </td><td> 0.00 </td><td> 0.00 </td><td> 0.55 </td><td> 3.62 </td><td> 0.000</td><td>⋯     </td><td> 0.00 </td><td> 0.41 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.48 </td><td> 0.200</td><td> 0.00 </td><td> 0.00 </td><td> 0.810</td></tr>
	<tr><th scope=row>uniprot:Q9Y6Z7</th><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td><td>⋯     </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td><td> 0.00 </td><td> 0.00 </td><td> 0.000</td></tr>
</tbody>
</table>



#### Filter only heart tissues columns

We got just two tissues to work with, one taken from a diseased heart and the second from other diseased heart but after a infarct.


```R
hearts <- grepl("heart", colnames(small_dts))
small_dts <- small_dts[hearts]
head(small_dts)
```


<table>
<thead><tr><th></th><th scope=col>heart..adult..diseased.post.infarction..donor1.CNhs11757.10050.101G5</th><th scope=col>heart..adult..diseased..donor1.CNhs11758.10051.101G6</th></tr></thead>
<tbody>
	<tr><th scope=row>uniprot:A0A183</th><td>0.2</td><td>0  </td></tr>
	<tr><th scope=row>uniprot:A0A577,uniprot:Q5ZGI7</th><td>0.0</td><td>0  </td></tr>
	<tr><th scope=row>uniprot:A0A582</th><td>0.0</td><td>0  </td></tr>
	<tr><th scope=row>uniprot:A0A583</th><td>0.1</td><td>0  </td></tr>
	<tr><th scope=row>uniprot:A0A597</th><td>0.0</td><td>0  </td></tr>
	<tr><th scope=row>uniprot:A0A599</th><td>0.0</td><td>0  </td></tr>
</tbody>
</table>




```R
colnames(small_dts) <- c("post_infarct", "pre_infarct")
# transform(small_dts, post_infarct = as.numeric(post_infarct), pre_infarct = as.numeric(pre_infarct))

```

### Compute MA

MA is a measure that consists on plotting the sum of the gen expressions of two observations on the difference of them.
Thereupon, 
- A = observ1 + observ2
- M = observ1 - observ2

The more similar the values are, the closer to the x axis they will be.
We use scale to normalize the data. It will calculate the mean and standard deviation (sd) of the entire vector, then "scale" each element by those values by subtracting the mean and dividing by the sd.


```R
A <- scale(rowSums(small_dts))
M <- scale(rowSums(transform(small_dts, pre_infarct = -(pre_infarct))))
head(M)
```


<table>
<tbody>
	<tr><th scope=row>uniprot:A0A183</th><td>0.03919355</td></tr>
	<tr><th scope=row>uniprot:A0A577,uniprot:Q5ZGI7</th><td>0.03746721</td></tr>
	<tr><th scope=row>uniprot:A0A582</th><td>0.03746721</td></tr>
	<tr><th scope=row>uniprot:A0A583</th><td>0.03833038</td></tr>
	<tr><th scope=row>uniprot:A0A597</th><td>0.03746721</td></tr>
	<tr><th scope=row>uniprot:A0A599</th><td>0.03746721</td></tr>
</tbody>
</table>



#### Plot results

We can observe that there are a few outliers.


```R
library(ggplot2)

sct <- ggplot(small_dts, aes(y=M, x=A)) + geom_point()
sct
```




![png](output_20_1.png)


We zoom on the two outliers on the lower right corner and show their uniprot_ids.


```R
sct + coord_cartesian(xlim = c(80,90),ylim= c(-75,-105)) +  geom_text(aes(label=rownames(small_dts)), size = 5)
```




![png](output_22_1.png)


And do the same on the center of the graphic where the are some values that differ, less than the two above, from the others.


```R
sct + coord_cartesian(xlim = c(25,60),ylim= c(-25,-50)) +  geom_text(aes(label=rownames(small_dts)), size = 3)
```




![png](output_24_1.png)


Unless we use an interactive plot, this does not seem to be a good method to visualize outliers.

### Compute K-Means


```R
# small_dts["uniprot:F5H4L2,uniprot:P12883","heart..adult..diseased.post.infarction..donor1.CNhs11757.10050.101G5"]
# uniprot:D3YTC6,uniprot:F6SZA2
```


```R
transp <- t(small_dts)
```

#### Matrix of distances.


```R
require(graphics)
mink_dist <- dist(transp, method='minkowski') # p,q >= 1
```

#### K-Means for 2 clusters.

We chose 2 clusters because we want to divide gens in two groups: One that affects on infarcts, on the other that does not.


```R
nhc_kmeans <- kmeans(small_dts,2)
```

Looking at the results, it seems that there are some genes that modify their behaviour once an infarct affects the patient.


```R
cat("Size of clusters")
data.frame(nhc_kmeans$size, row.names = c('Cluster 1', 'Cluster 2'))
```

    Size of clusters


<table>
<thead><tr><th></th><th scope=col>nhc_kmeans.size</th></tr></thead>
<tbody>
	<tr><th scope=row>Cluster 1</th><td>   15</td></tr>
	<tr><th scope=row>Cluster 2</th><td>31437</td></tr>
</tbody>
</table>



Now, we take the ids of the proteins within their cluster assignation and filter those that are in the cluster 1.


```R
htmap <- t(as.data.frame(nhc_kmeans$cluster))
fil <- htmap[,htmap[1,]==1, drop=FALSE]
colnames(fil)
```


<ol class=list-inline>
	<li>'uniprot:B0QYF7,uniprot:B4DUI1,uniprot:P02144'</li>
	<li>'uniprot:B4DL87,uniprot:F8WE04'</li>
	<li>'uniprot:D3YTC6,uniprot:F6SZA2'</li>
	<li>'uniprot:E7ET18,uniprot:A6NKB1,uniprot:E7EQE6,uniprot:Q8WZ42,uniprot:E9PPD3'</li>
	<li>'uniprot:F5H4L2,uniprot:P12883'</li>
	<li>'uniprot:F8W7I2,uniprot:P41222,uniprot:Q5SQ11'</li>
	<li>'uniprot:O14558,uniprot:E9PCI5'</li>
	<li>'uniprot:P05413'</li>
	<li>'uniprot:P10916,uniprot:G3V1V8'</li>
	<li>'uniprot:P17661,uniprot:Q53SB5'</li>
	<li>'uniprot:P63316,uniprot:Q6FH91'</li>
	<li>'uniprot:Q16654,uniprot:B3KU25,uniprot:A4D1H4'</li>
	<li>'uniprot:Q5T8M8,uniprot:Q5T8M7,uniprot:A6NL76,uniprot:P68133'</li>
	<li>'uniprot:Q9BUF6,uniprot:P45379,uniprot:F8WAF6,uniprot:Q7Z554'</li>
	<li>'uniprot:Q9Y427,uniprot:Q6ZN40,uniprot:P09493'</li>
</ol>



Now, we can take a look on the original values and see directly their variations.


```R
only_sp <- dts[dts$uniprot_id %in% colnames(fil),]
hearts <- grepl("heart", colnames(only_sp))
only_sp <- only_sp[hearts]
```


```R
rownames(only_sp) <- colnames(fil)
head(only_sp)
```


<table>
<thead><tr><th></th><th scope=col>heart..adult..diseased.post.infarction..donor1.CNhs11757.10050.101G5</th><th scope=col>heart..adult..diseased..donor1.CNhs11758.10051.101G6</th></tr></thead>
<tbody>
	<tr><th scope=row>uniprot:B0QYF7,uniprot:B4DUI1,uniprot:P02144</th><td>2664.52 </td><td> 6512.22</td></tr>
	<tr><th scope=row>uniprot:B4DL87,uniprot:F8WE04</th><td>7359.37 </td><td>19336.63</td></tr>
	<tr><th scope=row>uniprot:D3YTC6,uniprot:F6SZA2</th><td>2999.18 </td><td> 5711.28</td></tr>
	<tr><th scope=row>uniprot:E7ET18,uniprot:A6NKB1,uniprot:E7EQE6,uniprot:Q8WZ42,uniprot:E9PPD3</th><td>2290.52 </td><td> 5024.95</td></tr>
	<tr><th scope=row>uniprot:F5H4L2,uniprot:P12883</th><td>9245.94 </td><td>18377.74</td></tr>
	<tr><th scope=row>uniprot:F8W7I2,uniprot:P41222,uniprot:Q5SQ11</th><td>3167.91 </td><td> 7903.77</td></tr>
</tbody>
</table>



### Improving accuracy

Thanks to the first visualization on the MA plot, we saw that there were 2 points that differed more than the others.
We can then cluster over the previous cluster and achieve the same result.


```R
ds <- dist(only_sp)
disease_km <- kmeans(ds, 2)
```


```R
library(ggfortify)
set.seed(1)
clt <- kmeans(only_sp, 2)
autoplot(clt, data = only_sp, frame=T, frame.type = "convex")
```




![png](output_42_1.png)


We get the names of the outliers within the cluster 2


```R
htmap <- t(as.data.frame(clt$cluster))
fil <- htmap[,htmap[1,]==2, drop=FALSE]
colnames(fil)
```


<ol class=list-inline>
	<li>'uniprot:B4DL87,uniprot:F8WE04'</li>
	<li>'uniprot:F5H4L2,uniprot:P12883'</li>
</ol>



### Getting information about the proteins.

Now that we have identified two proteins that behave different, we just have to search about them to see what we have found.
The first pair, B4DL87 and F8WE04 are similar proteins, so this might me the reason of measure them together.
For the second pair, there is not information about F5H4L2, so it is not possible to confirm the same fact.

The role of HSP27 is to inhibit apoptosis whereas Myosin is involved on muscle contraction.

- B4DL87: http://www.uniprot.org/uniprot/B4DL87
- F8WE04 (HSP27): http://www.uniprot.org/uniprot/F8WE04

 - +info: https://en.wikipedia.org/wiki/Hsp27#Clinical_significance
    
- F5H4L2: http://www.uniprot.org/uniprot/F5H4L2

- P12883 (Myosin): http://www.uniprot.org/uniprot/P12883 
 - +info: https://en.wikipedia.org/wiki/MYH7


