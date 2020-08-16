1.numeric stability

solved by softmax function

```
v = exp(v - max(v))/sum(exp(v - max(v)))
```



2.multimodality issue

MCMC with informative prior不受该问题影响。

EM可能会受影响，当初始值选的很不好时，可能导致assign到某一个component数量为0：

进而导致w[k]=0，从而使得数值上为-inf.

```
log(w[k]) + xxx
```



3.BIC

trade-off of goodness of fit and complexity

