# Miller-Rabin Primality Test Experiments

This project implements the Miller-Rabin primality test and includes experiments to find composite numbers that pass the test incorrectly (false positives). The goal is to empirically study the error rate and compare it to the theoretical maximum of 1/4.

## Files

- `miller_rabin.py`: Core implementation of the Miller-Rabin primality test
- `experiment.py`: Script to run experiments, find false positives, and generate visualizations
- `requirements.txt`: Project dependencies
- `README.md`: This file

## Usage

First, install the required dependencies:

```bash
pip install -r requirements.txt
```

To run the experiments:

```bash
python experiment.py
```

This will:
1. Test different ranges of numbers for false positives
2. Include specific ranges containing known Carmichael numbers
3. Try different numbers of rounds (k) to see how it affects accuracy
4. Calculate error rates and compare with the theoretical 1/4 bound
5. Generate visualization plots:
   - `error_rate_vs_k.png`: Shows how error rate changes with number of rounds
   - `false_positives_vs_k.png`: Shows how false positives decrease with more rounds
   - `carmichael_success_rate.png`: Shows how well Carmichael numbers fool the test
6. Print a comprehensive summary including:
   - Overall results for each k value
   - Analysis of Carmichael numbers
   - Most successful false positives found

## About the Implementation

The Miller-Rabin test is a probabilistic primality test. For any composite number n, the probability that the test incorrectly identifies it as prime is at most 1/4 for each round. The probability decreases exponentially with more rounds.

Some interesting numbers to look for:
- 561 (smallest Carmichael number)
- 1105
- 1729 (Carmichael number)
- 2465
- 2821
- 6601

These numbers are known to be composite but can sometimes fool the Miller-Rabin test.

## Visualization Details

The experiment generates three main plots:

1. **Error Rate vs k**: Shows how the overall error rate decreases as we increase the number of rounds (k). Includes a reference line at 1/4 for comparison with the theoretical maximum.

2. **False Positives Count vs k**: Displays the total number of false positives found for each value of k on a logarithmic scale, demonstrating the exponential decrease in false positives with additional rounds.

3. **Carmichael Numbers Success Rate**: A bar chart showing how often specific Carmichael numbers successfully fool the test for different values of k, helping identify which numbers are particularly good at defeating the test.

## Requirements

- Python 3.6 or higher
- matplotlib >= 3.5.0
- numpy >= 1.21.0
- seaborn >= 0.11.0 