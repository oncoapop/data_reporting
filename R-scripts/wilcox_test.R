# Define data
# siRNA FACS from Takeda 72 hr
sg1<-c(7,9,12,17,19)
g1<-c(66,64,51,52,39)
s<-c(24,25,26,22,27)
g2m<-c(11,11,24,26,34)
A<-factor(c("siNT","siNT","sieIF4A3","sieIF4A3","sieIF4A3"))

# Define data# independent 2-group Mann-Whitney U Test 
wilcox.test(y~A) 
# where y is numeric and A is A binary factor

# independent 2-group t-test
t.test(sg1~A) # where y is numeric and x is a binary factor
t.test(g1~A) # where y is numeric and x is a binary factor
t.test(s~A) # where y is numeric and x is a binary factor
t.test(g2m~A) # where y is numeric and x is a binary factor

# Kruskal Wallis Test One Way Anova by Ranks 
subG1<-kruskal.test(sg1~A) # where y1 is numeric and A is a factor
G1G0<-kruskal.test(g1~A) # where y1 is numeric and A is a factor
S<-kruskal.test(s~A) # where y1 is numeric and A is a factor
G2M<-kruskal.test(g2m~A) # where y1 is numeric and A is a factor

subG1
G1G0
S
G2M


# independent 2-group Mann-Whitney U Test
y<-c(7,9)
x<-c(12,17,19)

wilcox.test(y,x) # where y and x are numeric

# independent 2-group t-test
t.test(x,y) # where y1 and y2 are numeric

