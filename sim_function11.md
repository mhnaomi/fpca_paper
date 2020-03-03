This file demonstrates the process of proposed hypothesis testing on simulated data.

The process requires packages 'MASS', 'nlme', 'lme4', and 'lmerTest'.

```
my_packages <- c("MASS", "nlme", "lme4", "lmerTest")
suppressMessages(lapply(my_packages, library, character.only = TRUE))
```

The functional data is simulated from the true model 
<img src="http://bit.ly/2PKvf2L" align="center" border="0" alt="Y_{gljv}(t)=\mu(t)+\tau_g(t)+U_{gl}(t) + W_{glj}(t) + \epsilon_{gljv}(t)" width="394" height="21" />,
where larger <img src="http://bit.ly/32N57tB" align="center" border="0" alt="b_\tau" width="22" height="18" /> corresponds to a larger difference between the two groups and 0 corresponds to no difference. 
