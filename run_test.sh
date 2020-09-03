# Some tests

cd main

# 1.  Quantifying SVs
num_sim=10000
/usr/bin/time python abc.py -N $num_sim > std_test1 &

# The plot will be available in the following directory: ../output/sumstats/sv_count.png



# 2.  Probability of Chromothripsis under Biased Chromosome Selection
num_sim=10000
odir=../output/sumstats/prob_trends/selection

for i in "random" "biased" "fixed"
do
  echo $i
  bias=$i
  prefix="prob_${bias}"
  /usr/bin/time python abc.py -a check_chromothripsis -N $num_sim -b $bias -o $odir -p $prefix > std_$prefix &
done

python chromothProb_selection.py




# 3.  Probability of Chromothripsis under Variable G1 repairs
num_sim=10000
odir=../output/sumstats/prob_trends/unconnected_G1

for i in 0 5 10 50 100
do
  echo $i
  lambda=$i
  prefix="lambda_${lambda}"
  /usr/bin/time python abc.py -a check_chromothripsis -N $num_sim -l $lambda -o $odir -p $prefix > std_$prefix &
done

python chromothProb_g1.py



# 4.  Probability of Chromothripsis under Variable nDSBs
num_sim=10000
odir=../output/sumstats/prob_trends/nDSBs

for i in 0 5 10 50 100
do
 echo $i
 minDSB=$i
 maxDSB=$i
 prefix="nDSB_${i}"
 /usr/bin/time python abc.py -a check_chromothripsis -N $num_sim -m $minDSB -M $maxDSB -o $odir -p $prefix > std_$prefix &
done

python chromothProb_nDSBs.py



# 5.  Probability of Chromothripsis under nCell Cycles
num_sim=10000
odir=../output/sumstats/prob_trends/nCycles/

for i in 1 2
do
 echo $i
 nCycles=$i
 prefix="${i}cycle"
 /usr/bin/time python abc.py -a check_chromothripsis -N $num_sim -n $nCycles -o $odir -p $prefix > std_$prefix &
done

python chromothProb_cycles.py
