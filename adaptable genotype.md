## GENOTYPE CALLING ASSIGNMENT - rs1126809 SNP

#### This python script is designed to create a VCF like file based on corresponding base calls and base qualities. The user must first input their data and then the script will run, creating tables for each step in the process before creating a VCF file like table out of the generated data. These steps are fully automated.


```python
#sequencer data input by user
BaseCalls = input('Enter base calls').upper()
BaseQualities = input('Enter base qualities').upper()

#reference and alt allele input by user
refallele = input('Enter reference allele').upper()
altallele = input('Enter alternate allele').upper()

#allele frequencies from dbSNP input by user
RefFreq = float(input('Enter frequency for reference allele'))
AltFreq = float(input('Enter frequency for alternative allele'))

#other information needed for the creation of teh VCF file like table input by user 
chromosome = input('Please enter on which chromosome this SNP lies:')
position = input('Please enter the position of the SNP on this chromosome:')
SNPID = input('Please enter the SNP ID:')
altallelefreq = 'AF:'+ str(AltFreq)
```

    Enter base calls CCCT
    Enter base qualities ?5+?
    Enter reference allele C
    Enter alternate allele T
    Enter frequency for reference allele 0.5
    Enter frequency for alternative allele 0.5
    Please enter on which chromosome this SNP lies: 1
    Please enter the position of the SNP on this chromosome: 1
    Please enter the SNP ID: 1



```python
#parameters for using genotype frequencies or flat priors
flatpriors = 0.33
refrefFreq = RefFreq*RefFreq
refaltFreq = 2*RefFreq*AltFreq
altaltFreq = AltFreq*AltFreq
```


```python
#convert strings to lists for indexing
def Convert(string):
    list1 = []
    list1[:0] = string
    return list1

BaseCalls = Convert(BaseCalls)
BaseQualities = Convert(BaseQualities)
```


```python
#install module tabulate to make tables
! pip3 install tabulate
from tabulate import tabulate
```

    Collecting tabulate
      Using cached https://files.pythonhosted.org/packages/92/4e/e5a13fdb3e6f81ce11893523ff289870c87c8f1f289a7369fb0e9840c3bb/tabulate-0.8.10-py3-none-any.whl
    Installing collected packages: tabulate
    Successfully installed tabulate-0.8.10



```python
#calculate Q scores from ascii given in basequalities and calculate error and accuracy

#set up data needed for table
col_names = ['Base', 'Q/Phred Score', 'Error', 'Accuracy']
Qerror = []

#cycle through each base call and calculate phred score error and accuracy 
for x in range(len(BaseQualities)):
    phredscoreconv = ord(BaseQualities[x])-33
    error = 10**(phredscoreconv/-10)
    accuracy = 1-error
    Qerror += [[BaseCalls[x], phredscoreconv, error, accuracy]]

print(tabulate(Qerror, headers=col_names, tablefmt="pretty"))
```

    +------+---------------+-------+----------+
    | Base | Q/Phred Score | Error | Accuracy |
    +------+---------------+-------+----------+
    |  C   |      30       | 0.001 |  0.999   |
    |  C   |      20       | 0.01  |   0.99   |
    |  C   |      10       |  0.1  |   0.9    |
    |  T   |      30       | 0.001 |  0.999   |
    +------+---------------+-------+----------+



```python
#create table for caluclating probability of base given the genotype 

#set up variables needed to create table 
baseprobscore = []
col_names2 = ['Base Probability', '0/0 Score', '0/1 Score', '1/1 Score']

#proability of data given genotype requires multiplication of all columns of each genotype
#set score to one and multiply this by each value by each entered value in the column
#this will calcualte probaility of base given genotype and probaility of data given genotype at the same time
refrefscoretotal = 1
refaltscoretotal = 1
altaltscoretotal = 1

#cycle trhough each base call and print different calcualted numbers based on how many matches of the base there were in the genotype
for x in range(len(BaseCalls)):
    
    #set up possible genotypes
    refref = refallele+refallele
    refalt = refallele+altallele
    altalt = altallele+altallele
    
    #define phredscores error and accuracy
    phredscoreconv = ord(BaseQualities[x])-33
    error = 10**(phredscoreconv/-10)
    accuracy = 1-error
    
    #bases match how many alleles formulae
    twomatch = accuracy 
    onematch = ((error/3)+accuracy)/2
    nomatch = error/3
        
    #make value entered into table for each base call dependent on number of matches of the base found in the genotype string
    if refref.count(BaseCalls[x])==0:
        refrefscore = nomatch
        refrefscoretotal *= refrefscore
    if refref.count(BaseCalls[x]) == 1:
        refrefscore = onematch
        refrefscoretotal *= refrefscore
    if refref.count(BaseCalls[x]) == 2:
        refrefscore = twomatch
        refrefscoretotal *= refrefscore
    
    if refalt.count(BaseCalls[x])==0:
        refaltscore = nomatch
        refaltscoretotal *= refaltscore
    if refalt.count(BaseCalls[x]) == 1:
        refaltscore = onematch
        refaltscoretotal *= refaltscore
    if refalt.count(BaseCalls[x]) == 2:
        refaltscore = twomatch
        refaltscoretotal *= refaltscore

    if altalt.count(BaseCalls[x]) == 0:
        altaltscore = nomatch
        altaltscoretotal *= altaltscore
    if altalt.count(BaseCalls[x]) == 1:
        altaltscore = onematch
        altaltscoretotal *= altaltscore
    if altalt.count(BaseCalls[x]) == 2:
        altaltscore = twomatch
        altaltscoretotal *= altaltscore
    
    #append a new row into the table for every base call with its corresponding value for GG GA and AA genotypes
    baseprobscore+=[[BaseCalls[x], refrefscore, refaltscore, altaltscore]]
    
#print table of base probability given the genotype
print('Base probability scores: \n')
print(tabulate(baseprobscore, headers=col_names2, tablefmt="pretty"))
```

    Base probability scores: 
    
    +------------------+-----------------------+---------------------+-----------------------+
    | Base Probability |       0/0 Score       |      0/1 Score      |       1/1 Score       |
    +------------------+-----------------------+---------------------+-----------------------+
    |        C         |         0.999         | 0.49966666666666665 | 0.0003333333333333333 |
    |        C         |         0.99          | 0.49666666666666665 | 0.0033333333333333335 |
    |        C         |          0.9          | 0.4666666666666667  |  0.03333333333333333  |
    |        T         | 0.0003333333333333333 | 0.49966666666666665 |         0.999         |
    +------------------+-----------------------+---------------------+-----------------------+



```python
#print table of P(DATA|GENOTYPE)
dataprobscore = [['P(DATA|GENOTYPE)', refrefscoretotal, refaltscoretotal, altaltscoretotal]]
col_names3 = ['GENOTYPE', '0/0', '0/1', '1/1']
print('P(DATA|GENOTYPE): \n')
print(tabulate(dataprobscore, col_names3, tablefmt="pretty"))
```

    P(DATA|GENOTYPE): 
    
    +------------------+------------------------+--------------------+-----------------------+
    |     GENOTYPE     |          0/0           |        0/1         |          1/1          |
    +------------------+------------------------+--------------------+-----------------------+
    | P(DATA|GENOTYPE) | 0.00029670299999999994 | 0.0578672109382716 | 3.700000000000001e-08 |
    +------------------+------------------------+--------------------+-----------------------+



```python
#calculating genotype posterior probabilities
#P(genotype|data)*P(data)=P(DATA|GENOTYPE)*P(GENOTYPE)
#P(data)=sum of 3 values for above 

#USING FLAT PRIORS
#P(GENOTYPE|DATA)*P(GENOTYPE)
refrefposteriorflat = refrefscoretotal*flatpriors
refaltposteriorflat = refaltscoretotal*flatpriors
altaltposteriorflat = altaltscoretotal*flatpriors

#USING GENOTYPE FREQUENCY 
refrefposteriorGF = refrefscoretotal*refrefFreq
refaltposteriorGF = refaltscoretotal*refaltFreq
altaltposteriorGF = altaltscoretotal*altaltFreq

pdataflat = refrefposteriorflat + refaltposteriorflat + altaltposteriorflat
pdataGF = refrefposteriorGF + refaltposteriorGF + altaltposteriorGF

print('P(DATA|GENOTYPE)*P(GENOTYPE) Using flat priors: \n')
posteriorprobflat = [['P(GENOTYPE|DATA)*P(DATA)', refrefposteriorflat, refaltposteriorflat, altaltposteriorflat, pdataflat]]
col_names4 = ['GENOTYPE', '0/0', '0/1', '1/1', 'P(data)']
print(tabulate(posteriorprobflat, col_names4, tablefmt="pretty"), '\n\n')

print('P(DATA|GENOTYPE)*P(GENOTYPE) Using genotype frequencies: \n')
posteriorprobGF = [['P(GENOTYPE|DATA)*P(DATA)', refrefposteriorGF, refaltposteriorGF, altaltposteriorGF, pdataGF]]
print(tabulate(posteriorprobGF, col_names4, tablefmt="pretty"))
```

    P(DATA|GENOTYPE)*P(GENOTYPE) Using flat priors: 
    
    +--------------------------+-----------------------+---------------------+------------------------+---------------------+
    |         GENOTYPE         |          0/0          |         0/1         |          1/1           |       P(data)       |
    +--------------------------+-----------------------+---------------------+------------------------+---------------------+
    | P(GENOTYPE|DATA)*P(DATA) | 9.791198999999998e-05 | 0.01909617960962963 | 1.2210000000000003e-08 | 0.01919410380962963 |
    +--------------------------+-----------------------+---------------------+------------------------+---------------------+ 
    
    
    P(DATA|GENOTYPE)*P(GENOTYPE) Using genotype frequencies: 
    
    +--------------------------+-----------------------+--------------------+-----------------------+----------------------+
    |         GENOTYPE         |          0/0          |        0/1         |          1/1          |       P(data)        |
    +--------------------------+-----------------------+--------------------+-----------------------+----------------------+
    | P(GENOTYPE|DATA)*P(DATA) | 7.417574999999999e-05 | 0.0289336054691358 | 9.250000000000002e-09 | 0.029007790469135798 |
    +--------------------------+-----------------------+--------------------+-----------------------+----------------------+



```python
#PROVE ALL GENOTYPE POSTERIOR PROBABILITIES ADD UP TO ONE TO SOLVE FOR P(DATA)
#p(data) flat priors
datarefref = refrefposteriorflat/pdataflat
datarefalt = refaltposteriorflat/pdataflat
dataaltalt = altaltposteriorflat/pdataflat
sumpdata = datarefref + datarefalt + dataaltalt
col_names5 = ['GENOTYPE', '0/0', '0/1', '1/1']
genposdata = [['P(GENOTYPE|DATA)', datarefref, datarefalt, dataaltalt]]
print(tabulate(genposdata, col_names5, tablefmt="pretty"))
print('Sum of genotype posterior probabilities using flat posteriors:', sumpdata, '\n\n')

#p(data) genotype frequencies 
GFdatarefref = refrefposteriorGF/pdataGF
GFdatarefalt = refaltposteriorGF/pdataGF
GFdataaltalt = altaltposteriorGF/pdataGF
GFsumpdata = GFdatarefref + GFdatarefalt + GFdataaltalt
GFgenposdata = [['P(GENOTYPE|DATA)', GFdatarefref, GFdatarefalt, GFdataaltalt]]
print(tabulate(GFgenposdata, col_names5, tablefmt="pretty"))
print('Sum of genotype posterior probabilities using  genotype frequencies:', GFsumpdata)
```

    +------------------+----------------------+--------------------+-----------------------+
    |     GENOTYPE     |         0/0          |        0/1         |          1/1          |
    +------------------+----------------------+--------------------+-----------------------+
    | P(GENOTYPE|DATA) | 0.005101149341021996 | 0.9948982145261258 | 6.361328521040028e-07 |
    +------------------+----------------------+--------------------+-----------------------+
    Sum of genotype posterior probabilities using flat posteriors: 0.9999999999999999 
    
    
    +------------------+-----------------------+--------------------+-----------------------+
    |     GENOTYPE     |          0/0          |        0/1         |          1/1          |
    +------------------+-----------------------+--------------------+-----------------------+
    | P(GENOTYPE|DATA) | 0.0025570975520842498 | 0.9974425835680615 | 3.188798543564348e-07 |
    +------------------+-----------------------+--------------------+-----------------------+
    Sum of genotype posterior probabilities using  genotype frequencies: 1.0



```python
#raw PLs
import math as m

#flat
refrefflatPL = round(-10*m.log(datarefref, 10), 0)
refaltflatPL = round(-10*m.log(datarefalt, 10), 0)
altaltflatPL = round(-10*m.log(dataaltalt, 10), 0)
col_names6 = ['GENOTYPE', '0/0', '0/1', '1/1']
rawplflat = [['RAW PL(GENOTYPE)', refrefflatPL, refaltflatPL, altaltflatPL]]
print('Raw PL values for flat posteriors:')
print(tabulate(rawplflat, col_names6, tablefmt="pretty"), '\n\n')

#GF
refrefPLgf = round(-10*m.log(GFdatarefref, 10), 0)
refaltPLgf = round(-10*m.log(GFdatarefalt, 10), 0)
altaltPLgf = round(-10*m.log(GFdataaltalt, 10), 0)
rawplgf = [['RAW PL(GENOTYPE)', refrefPLgf, refaltPLgf, altaltPLgf]]
print('Raw PL values for genotype frequencies')
print(tabulate(rawplgf, col_names6, tablefmt="pretty"), '\n\n')

#normalise PL values by taking lowest value away from each raw PL
minflat = min(refrefflatPL, refaltflatPL, altaltflatPL)
normalisedrefrefflatPL = round(-10*m.log(datarefref, 10) - minflat, 0)
normalisedrefaltflatPL = round(-10*m.log(datarefalt, 10) - minflat, 0)
normalisedaltaltflatPL = round(-10*m.log(dataaltalt, 10) - minflat, 0)
normalisedplflat = [['RAW PL(GENOTYPE)', normalisedrefrefflatPL, normalisedrefaltflatPL, normalisedaltaltflatPL]]
print('Normalised PL values for flat posteriors:')
print(tabulate(normalisedplflat, col_names6, tablefmt="pretty"), '\n\n')

minGF = min(refrefPLgf, refaltPLgf, altaltPLgf)
normalisedrefrefGFPL = round(-10*m.log(GFdatarefref, 10) - minGF, 0)
normalisedrefaltGFPL = round(-10*m.log(GFdatarefalt, 10) - minGF, 0)
normalisedaltaltGFPL = round(-10*m.log(GFdataaltalt, 10) - minGF, 0)
normalisedplgf = [['RAW PL(GENOTYPE)', normalisedrefrefGFPL, normalisedrefaltGFPL, normalisedaltaltGFPL]]
print('Normalised PL values for genotype frequencies:')
print(tabulate(normalisedplgf, col_names6, tablefmt="pretty"), '\n\n')

```

    Raw PL values for flat posteriors:
    +------------------+------+-----+------+
    |     GENOTYPE     | 0/0  | 0/1 | 1/1  |
    +------------------+------+-----+------+
    | RAW PL(GENOTYPE) | 23.0 | 0.0 | 62.0 |
    +------------------+------+-----+------+ 
    
    
    Raw PL values for genotype frequencies
    +------------------+------+-----+------+
    |     GENOTYPE     | 0/0  | 0/1 | 1/1  |
    +------------------+------+-----+------+
    | RAW PL(GENOTYPE) | 26.0 | 0.0 | 65.0 |
    +------------------+------+-----+------+ 
    
    
    Normalised PL values for flat posteriors:
    +------------------+------+-----+------+
    |     GENOTYPE     | 0/0  | 0/1 | 1/1  |
    +------------------+------+-----+------+
    | RAW PL(GENOTYPE) | 23.0 | 0.0 | 62.0 |
    +------------------+------+-----+------+ 
    
    
    Normalised PL values for genotype frequencies:
    +------------------+------+-----+------+
    |     GENOTYPE     | 0/0  | 0/1 | 1/1  |
    +------------------+------+-----+------+
    | RAW PL(GENOTYPE) | 26.0 | 0.0 | 65.0 |
    +------------------+------+-----+------+ 
    
    



```python
#genotype quality 
import statistics

#flat
flatGQlist = [refrefflatPL, refaltflatPL, altaltflatPL]
flatGQ = statistics.median(flatGQlist)
print('For flat posteriors the genotype quality is:', flatGQ)

#genotype frequencies
GFGQlist = [refrefPLgf, refaltPLgf, altaltPLgf]
GFGQ = statistics.median(GFGQlist)
print('For genotype frequencies the genotype quality is:', GFGQ)
```

    For flat posteriors the genotype quality is: 23.0
    For genotype frequencies the genotype quality is: 26.0



```python
#creating a variable for full stramore man sample information 
#tabulate cannot accomodate more than one varibale in a cell at a time and so a varibale must be created outside the table containing all information

#find lowest PL score and which genotype it belongs to
#flat
FLATrefrefplmin = normalisedrefrefflatPL
FLATrefaltplmin = normalisedrefaltflatPL
FLATaltaltplmin = normalisedaltaltflatPL
minPLflat = {FLATrefrefplmin:"0/0",FLATrefaltplmin:"0/1",FLATaltaltplmin:"1/1"}
minPLflatvar = minPLflat.get(min(minPLflat))
minPLflatnumber = min(minPLflat)

#GF
GFrefrefplmin = normalisedrefrefGFPL
GFrefaltplmin = normalisedrefaltGFPL
GFaltaltplmin = normalisedaltaltGFPL
minPLGF = {GFrefrefplmin:"0/0",GFrefaltplmin:"0/1",GFaltaltplmin:"1/1"}
minPLGFvar = minPLGF.get(min(minPLGF))
minPLGFnumber = min(minPLGF)

#flat sramore variable 
flatsramorevar = minPLflatvar + ':' + str(int(flatGQ)) + ':' + str(len(BaseCalls)) + ':' + str(int(minPLflatnumber))

#genotype frequency sramore variable 
GFsramorevar = minPLGFvar + ':' + str(int(GFGQ)) + ':' + str(len(BaseCalls)) + ':' + str(int(minPLGFnumber))
```


```python
#VCF LIKE FILE CREATION 

VCFHEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Sramore_man']

VCFFLAT = [[chromosome, position, SNPID, refallele, altallele, '.', '.', '.', 'GT:GQ:DP:PL', flatsramorevar]]

print('VCF FILE USING FLAT PROBABILITIES:')
print(tabulate(VCFFLAT, VCFHEADER, tablefmt="pretty"), '\n\n')

VCFGF = [[chromosome, position, SNPID, refallele, altallele, '.', '.', altallelefreq, 'GT:GQ:DP:PL', GFsramorevar]]

print('VCF FILE USING GENOTYPE FREQUENCIES:')
print(tabulate(VCFGF, VCFHEADER, tablefmt="pretty"))
```

    VCF FILE USING FLAT PROBABILITIES:
    +-------+-----+----+-----+-----+------+--------+------+-------------+-------------+
    | CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO |   FORMAT    | Sramore_man |
    +-------+-----+----+-----+-----+------+--------+------+-------------+-------------+
    |   1   |  1  | 1  |  C  |  T  |  .   |   .    |  .   | GT:GQ:DP:PL | 0/1:23:4:0  |
    +-------+-----+----+-----+-----+------+--------+------+-------------+-------------+ 
    
    
    VCF FILE USING GENOTYPE FREQUENCIES:
    +-------+-----+----+-----+-----+------+--------+--------+-------------+-------------+
    | CHROM | POS | ID | REF | ALT | QUAL | FILTER |  INFO  |   FORMAT    | Sramore_man |
    +-------+-----+----+-----+-----+------+--------+--------+-------------+-------------+
    |   1   |  1  | 1  |  C  |  T  |  .   |   .    | AF:0.5 | GT:GQ:DP:PL | 0/1:26:4:0  |
    +-------+-----+----+-----+-----+------+--------+--------+-------------+-------------+

