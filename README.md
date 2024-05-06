# pept-prot-infer
`pept-prot-infer` determines the original peptide identities and their relative abundances from single molecule sequencing datasets, particularly fluorosequencing data preprocessed by https://github.com/marcottelab/whatprot. This peptide distribution is then used in a second step to figure out the original protein identities and abundances.

## Why determining proteins or even just peptides from SM-sequencing data is tricky
In addition to the physico-chemical problems that can happen during SM-sequencing (Edman degradation failure, photobleaching, the fact that many amino acids are not labelable, erroneous or no attachment of reporter labels, FRET, and more), there is also the simple statistical problem that there are only very few reads per peptide/protein in a large pool of peptides/proteins, as opposed to sequencing techniques in which one uses a large quantity of purified protein. Each read unfortunately does not identify a unique peptide, but could be be a wide slew of peptides, some with bigger probabilites, others with smaller ones.

Therefore, calculating even just the abundance of even just the peptides is quite tricky, and deriving the protein distribution is even more difficult.
Instead, it is easier to compute P(Read<sub>i</sub>|Peptide<sub>j</sub>), which is the likelihood for a read i to be generated assuming that a peptide j was present. 
[whatprot](https://github.com/marcottelab/whatprot), for example, uses an HMM approach to do this. Even with P(Read<sub>i</sub>|Peptide<sub>j</sub>), it is difficult to calculate P(Peptide<sub>j</sub>) directly. One may think that it'd simply be the number of partial sequences divided by the number of reads, but that is only accurate if the number of reads is much larger than the number of possible partial sequences. In SM-sequencing, this is typically not the case. `pept-prot-infer` instead, in the peptide inference step, converges on P(Peptide<sub>j</sub>) by using an expectation-maximization algorithm which alternates between approximating P(Peptide<sub>j</sub>|Read<sub>i</sub>) and P(Peptide<sub>j</sub>). In the protein inference step, `pept-prot-infer` uses this P(Peptide<sub>j</sub>) distribution to calculate the protein distribution P(Protein<sub>k</sub>), using the same algorithm.

## Prerequisites
Anything that can run .ipynb or .py files. I used 
 JupyterLab or VSC with the extensions for Jupyter or Python.

## Usage tutorial
The folder `/example input` contains an exemplatory SM-sequencing dataset with 93197 reads from a run with 10 original proteins and 110 detected digest peptides.

### Step 0: input files
Each row in `1 peptide infer input p_xi_given_zj.csv` is one read, containing the likelihood scores P(Peptide<sub>i</sub>|Read<sub>j</sub>); i.e. the likelihood of each peptide i given this read j. For calculating P(Peptide<sub>i</sub>|Read<sub>j</sub>) from SM-sequencing reads, you can use [whatprot](https://github.com/marcottelab/whatprot).

|            | Peptide 0 | Peptide 1 | Peptide 2 | Peptide 3 | ... | Peptide 109 |
|-----------:|----------:|----------:|----------:|----------:|-----|-------------|
|     Read 0 | 3.76E+29  | 3.76E+29  | 6.81E+29  | 3.76E+29  | ... | 0           |
|     Read 1 | 1.60E+23  | 1.60E+23  | 1.05E+23  | 1.60E+23  | ... | 0           |
|     Read 2 | 1.44E+28  | 7.37E+27  | 9.26E+23  | 7.37E+27  | ... | 0           |
| ...        | ...       | ...       | ...       | ...       | ... | ...         |
| Read 93196 | 3.76E+18  | 7.34E+18  | 2.35E+15  | 9.94E+17  | ... | 0           |

### Step 1: Peptide inferral
This calculates the relative abundancies of the digest peptides in the assay, P(peptide<sub>j</sub>). The resultant data is used for step 2, the protein inferral. The script may complain about the missing files for step 2, just ignore that for now.

Run the `pept-prot-infer.ipynb` script. User is asked setting paramenters before running (see table below). 

|Setting|Meaning|Default|
|:-|:-|:-|
|EM_convergence_minimum|The lower, the stricter. Checks for difference between previous and current run|0.0001|
|EM_loopcounter_max|Maximum number of EM runs (per bootstrap run)|200|
|bootstrap_sampled_fraction|Fraction of subarray sampled for each bootstrap run. I haven't figured out what the ideal is, but this seems good.|0.8|
|n_bootstrap_runs|Number of bootstrap runs. More = better accuracy, up to some point.|200|
|CI_percent|Confidence interval|95%|

Result:
| Peptide ID | AVG [%] | ±STD [%] | -CI [%] | +CI [%] |
|-----------:|--------:|---------:|--------:|--------:|
|          0 |    3.11 |     0.12 |    2.87 |    3.30 |
|          1 |    4.66 |     0.09 |    4.53 |    4.79 |
|          2 |    0.58 |     0.05 |    0.50 |    0.65 |
|         ...|      ...|       ...|      ...|      ...|
|         109|     0.15|	  0.03| 	0.11|     0.20|

![peptide inferral step result visualised](https://raw.githubusercontent.com/SophiaBailingZhou/pept-prot-infer/tree/main/.github/Peptide_conc._110pept10prot.emf)

To calculate `2a protein infer input dyeseq in prot match count.csv` and `2b protein infer input p_xi_given_zj.csv`, please perform a virtual digest and virtual labelling of the organism's proteome (or use a smaller custom list of possible proteins, if you can narrow it down even more) using the proteases and chemical labels you used in your assay. P(peptide<sub>j</sub>|protein<sub>k</sub>) can be calculated by first matching the peptides found in the peptide inferrand then dividing the the number matches of a peptide to one digested protein by the sum of matches all peptides have to that protein. Save this in `2b protein infer input p_xi_given_zj.csv`. `2a protein infer input dyeseq in prot match count.csv` contains the sum of the count of digested peptides per protein.
Run the script again. It will now also give you give you protein distributions.

![protein inferral step result visualised](https://raw.githubusercontent.com/SophiaBailingZhou/pept-prot-infer/tree/main/.github/Protein_conc._110pept10prot.emf)

## Explanation for what this code does in detail
**Part 1: Calculation of peptide distribution P(Zj) using Matt's P(Xi|Zj) scores**
1. Imports P(Xi|Zj) table and asks user for settings
2. For each bootstrap run:
    1. Create a subtable created from randomly sampling rows in P(Xi|Zj), with replacement (bootstrapping)
    2. Calculate P(Zj) via EM using P(Xi|Zj) (subtable)            
       - This is a classic "Hen or egg" problem. There are two missing variables: P(Zj|Xi) and P(Zj). Calculating either requires the other.
       - For 1st calculation of P(Zj|Xi), P(Zj) is approximated by assuming equal distribution
       - Then alternates between updating P(Zj|Xi) and P(Zj) until convergence

3. Calculate average P(Zj) from all bootstrap runs, and calculate STD, 95% CI**

**Part 2a: Virtual digest of target proteome and match finding**

*This part is not implemented in this code, but I imagine it should be pretty easy to implement.*
1. Virtually digest target proteome into peptides using the same proteases
2. Count number of matches for all peptides inferred from part 1 against digested proteome (save in `dyeseq_in_prot_match_count.csv`)
3. Calculate P(Xi|Zj) (likelihood of peptide Xi given protein Zj) by dividing the match count of each peptide in a each digested protein by the total number of match counts of ALL peptides against that digested protein

**Part 2b: Calculation of protein distribution using calculated peptide distribution**
1. run bootstrap-EM, like in peptide inferring step
    - First calculation of P(Zj|Xi) assumes equal distribution of all P(Zj)
    - Uses peptide distribution P(Zj) values from (1) as P(Xi), which simplifies the Bayes-derived equation. However, the principle is the same

## Ideas for improvement

- Get rid of bootstrap while loop through broadcasting
- Get rid of checking for EM convergence every EM run... Instead its probably faster to just e.g. run 199 EM runs, then check for convergence between run 199 and 200. If convergence is still not met, suggest increasing runs to user or just auto-run another 199 runs


## Paper demonstrating capabilities
<mark>link to paper</mark>