"""
Experiment script to find composite numbers that fool the Miller-Rabin test,
calculate error rates, and visualize results.
"""

from miller_rabin import find_false_positives, is_prime_mr, is_definitely_prime
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Dict
import time
from itertools import combinations
from datetime import datetime

def write_to_log(message: str, file):
    """Helper function to write to both console and log file."""
    print(message)
    file.write(message + "\n")

def save_results(all_results: List[Dict], carmichael_results: Dict[int, Dict[int, float]]) -> None:
    """
    Save all experimental results to a detailed log file.
    Only show range results with non-zero error rates.
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"miller_rabin_results_{timestamp}.txt"
    
    with open(filename, "w") as f:
        # Write header
        write_to_log("=" * 80, f)
        write_to_log("MILLER-RABIN PRIMALITY TEST EXPERIMENT RESULTS", f)
        write_to_log(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", f)
        write_to_log("=" * 80 + "\n", f)
        
        # Write experiment parameters
        write_to_log("EXPERIMENT PARAMETERS", f)
        write_to_log("-" * 40, f)
        write_to_log(f"k values tested: {[result['k'] for result in all_results]}", f)
        write_to_log(f"Number of Carmichael numbers tested: {len(next(iter(carmichael_results.values())))}", f)
        write_to_log(f"Total ranges tested: {len(all_results[0]['range_data'])}", f)
        write_to_log("\n" + "=" * 80 + "\n", f)
        
        # Write overall results for each k
        write_to_log("OVERALL RESULTS BY NUMBER OF ROUNDS", f)
        write_to_log("-" * 40, f)
        for result in all_results:
            k = result['k']
            write_to_log(f"\nk = {k}:", f)
            write_to_log(f"  Total composites tested: {result['total_composites']:,}", f)
            
            # Only show false positives and error rate if non-zero
            if result['total_false_positives'] > 0:
                write_to_log(f"  False positives found: {result['total_false_positives']:,}", f)
                write_to_log(f"  Overall error rate: {result['overall_error_rate']:.8f}", f)
                
                # Only show ranges with false positives
                interesting_ranges = [r for r in result['range_data'] if r['false_positives']]
                if interesting_ranges:
                    write_to_log("\n  Ranges with false positives:", f)
                    for range_data in interesting_ranges:
                        start, end = range_data['range']
                        write_to_log(f"\n  Range {start:,} to {end:,}:", f)
                        write_to_log(f"    Composites tested: {range_data['composites']:,}", f)
                        write_to_log(f"    False positives found: {len(range_data['false_positives'])}", f)
                        write_to_log(f"    Error rate: {range_data['error_rate']:.8f}", f)
                        write_to_log(f"    Examples: {', '.join(str(n) for n in range_data['false_positives'][:5])}", f)
            else:
                write_to_log("  No false positives found", f)
        
        write_to_log("\n" + "=" * 80 + "\n", f)
        
        # Write Carmichael numbers analysis (only show if they fooled the test)
        write_to_log("CARMICHAEL NUMBERS ANALYSIS", f)
        write_to_log("-" * 40, f)
        carmichael_numbers = sorted(list(carmichael_results[list(carmichael_results.keys())[0]].keys()))
        
        # Table header
        header = "Number    " + "    ".join(f"k={k:<6}" for k in sorted(carmichael_results.keys()))
        write_to_log("\n" + header, f)
        write_to_log("-" * len(header), f)
        
        # Table content - only show numbers that fooled the test at least once
        for n in carmichael_numbers:
            success_rates = [carmichael_results[k][n] for k in sorted(carmichael_results.keys())]
            if any(rate > 0 for rate in success_rates):
                row = f"{n:<9} "
                row += "    ".join(f"{rate*100:6.2f}%" for rate in success_rates)
                write_to_log(row, f)
        
        write_to_log("\n" + "=" * 80 + "\n", f)
        
        # Write most successful false positives
        write_to_log("MOST SUCCESSFUL FALSE POSITIVES", f)
        write_to_log("-" * 40, f)
        for result in all_results:
            k = result['k']
            all_fps = []
            for range_data in result['range_data']:
                all_fps.extend(range_data['false_positives'])
            
            if all_fps:
                # Sort by value and get unique numbers
                unique_fps = sorted(set(all_fps))
                write_to_log(f"\nk = {k}:", f)
                write_to_log(f"  Found {len(unique_fps)} unique false positives", f)
                write_to_log("  Top 10 examples: " + ", ".join(str(n) for n in unique_fps[:10]), f)
        
        write_to_log("\n" + "=" * 80, f)
        write_to_log("\nExperiment completed successfully.", f)
        write_to_log("Visualization plots have been saved:", f)
        write_to_log("1. error_rate_vs_k.png - Error rate vs number of rounds", f)
        write_to_log("2. false_positives_vs_k.png - False positives vs number of rounds", f)
        write_to_log("3. carmichael_success_rate.png - Carmichael numbers success rates", f)
    
    print(f"\nDetailed results have been saved to: {filename}")

def generate_complex_ranges() -> List[Tuple[int, int]]:
    """
    Generate test ranges that are likely to produce false positives.
    Includes ranges around Carmichael numbers and their products.
    """
    # Extended list of Carmichael numbers
    carmichael = [
        561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341,
        41041, 46657, 52633, 62745, 63973, 75361, 101101, 115921, 126217,
        162401, 172081, 188461, 252601, 278545, 294409, 314821, 334153,
        340561, 399001, 410041, 449065, 488881, 512461
    ]
    
    ranges = []
    
    # Test around each Carmichael number with wider ranges
    for c in carmichael:
        # Test range around the Carmichael number
        ranges.append((max(2, c - 100), c + 100))
        
        # Test range around twice the Carmichael number
        ranges.append((2*c - 100, 2*c + 100))
        
        # Test range around three times the Carmichael number
        ranges.append((3*c - 100, 3*c + 100))
    
    # Test products of pairs of Carmichael numbers
    small_carmichael = carmichael[:8]  # Use first 8 for products
    for c1, c2 in combinations(small_carmichael, 2):
        product = c1 * c2
        ranges.append((product - 200, product + 200))
    
    # Add more general ranges with better distribution
    ranges.extend([
        (1000, 5000), 
        (9900, 10100),
        (29340, 29342),        # Around specific Carmichael number 29341
        (50000, 51000),  
        (99000, 100000), 
        (561*1105, 561*1105+1000),  # Range after product of first two Carmichael numbers
        (1000000, 1001000),    
        (10000000, 10001000),  
        (100000000, 100001000), 
    ])
    
    # Add ranges around powers of 2
    for power in range(10, 30):
        base = 2 ** power
        ranges.append((base - 100, base + 100))
    
    # Add ranges around powers of 10
    for power in range(3, 9):
        base = 10 ** power
        ranges.append((base - 100, base + 100))
    
    return ranges

def run_experiment(ranges: List[Tuple[int, int]], k: int = 20) -> Dict:
    """
    Run experiments over multiple ranges to find false positives and calculate error rates.
    Returns the collected data for visualization.
    """
    total_composites = 0
    total_false_positives = 0
    range_data = []
    
    print(f"\nRunning Miller-Rabin experiments with k={k} rounds:")
    print("-" * 60)
    
    # Run each range multiple times for better statistical significance
    num_runs = 1000
    
    for start, end in ranges:
        print(f"\nTesting range {start:,} to {end:,}")
        
        # Store results from multiple runs
        run_results = []
        for run in range(num_runs):
            # Find false positives in this range
            false_positives = find_false_positives(start, end, k)
            run_results.append(false_positives)
        
        # Combine results from all runs
        all_false_positives = set()
        for result in run_results:
            all_false_positives.update(result)
        
        # Count composites in this range
        composites = sum(1 for n in range(start, end + 1) 
                        if not is_definitely_prime(n))
        
        # Calculate error rate for this range
        error_rate = len(all_false_positives) / composites if composites > 0 else 0
        
        print(f"Composites tested: {composites:,}")
        print(f"False positives found across {num_runs} runs: {len(all_false_positives)}")
        if all_false_positives:
            print("Examples of false positives:", 
                  ", ".join(str(n) for n in sorted(all_false_positives)[:5]))
        print(f"Error rate: {error_rate:.8f}")
        
        range_data.append({
            'range': (start, end),
            'composites': composites,
            'false_positives': list(all_false_positives),
            'error_rate': error_rate,
            'num_runs': num_runs
        })
        
        total_composites += composites
        total_false_positives += len(all_false_positives)
    
    overall_error_rate = (total_false_positives / total_composites 
                         if total_composites > 0 else 0)
    
    print("\nOverall Results:")
    print("-" * 60)
    print(f"Total composites tested: {total_composites:,}")
    print(f"Total false positives found: {total_false_positives}")
    print(f"Overall error rate: {overall_error_rate:.8f}")
    print(f"Theoretical maximum error rate: 0.25 (1/4)")
    
    return {
        'k': k,
        'total_composites': total_composites,
        'total_false_positives': total_false_positives,
        'overall_error_rate': overall_error_rate,
        'range_data': range_data,
        'num_runs': num_runs
    }

def test_carmichael_numbers(k: int) -> Dict[int, float]:
    """
    Test how well Carmichael numbers fool the Miller-Rabin test.
    Returns a dictionary mapping each Carmichael number to its success rate.
    """
    # Extended list of Carmichael numbers
    carmichael_numbers = [
        561, 1105, 1729, 2465, 2821, 6601, 8911, 10585, 15841, 29341,
        41041, 46657, 52633, 62745, 63973, 75361, 101101, 115921, 126217,
        162401, 172081, 188461, 252601, 278545, 294409, 314821, 334153,
        340561, 399001, 410041, 449065, 488881, 512461
    ]
    results = {}
    
    print(f"\nTesting Carmichael numbers with k={k}:")
    print("-" * 60)
    
    for n in carmichael_numbers:
        # Increase trials for more accurate results
        trials = 1000
        fooled_count = sum(1 for _ in range(trials) 
                          if is_prime_mr(n, k))
        fooled_rate = fooled_count / trials
        results[n] = fooled_rate
        
        print(f"Number {n}: fooled the test {fooled_rate*100:.1f}% of the time")
    
    return results

def create_visualizations(all_results: List[Dict], carmichael_results: Dict[int, Dict[int, float]]) -> None:
    """
    Create and save visualization plots.
    """
    # Set up plot style
    plt.rcParams['figure.figsize'] = [10, 6]
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.alpha'] = 0.3
    
    # Plot 1: Error Rate vs k
    plt.figure()
    plt.plot([result['k'] for result in all_results],
             [result['overall_error_rate'] for result in all_results],
             'bo-', label='Empirical Error Rate')
    plt.axhline(y=0.25, color='r', linestyle='--', label='Theoretical Maximum (1/4)')
    plt.xlabel('Number of Rounds (k)')
    plt.ylabel('Error Rate')
    plt.title('Miller-Rabin Test Error Rate vs Number of Rounds')
    plt.legend()
    plt.savefig('error_rate_vs_k.png')
    plt.close()

    # Plot 2: False Positives Count vs k
    plt.figure()
    plt.plot([result['k'] for result in all_results],
             [result['total_false_positives'] for result in all_results],
             'go-')
    plt.xlabel('Number of Rounds (k)')
    plt.ylabel('Number of False Positives')
    plt.title('Number of False Positives vs Number of Rounds')
    plt.yscale('log')
    plt.savefig('false_positives_vs_k.png')
    plt.close()

    # Plot 3: Carmichael Numbers Success Rate
    plt.figure(figsize=(12, 6))
    carmichael_numbers = sorted(list(carmichael_results[list(carmichael_results.keys())[0]].keys()))
    x = np.arange(len(carmichael_numbers))
    width = 0.15
    
    for i, k in enumerate(sorted(carmichael_results.keys())):
        rates = [carmichael_results[k][n] for n in carmichael_numbers]
        plt.bar(x + i*width, rates, width, label=f'k={k}')
    
    plt.xlabel('Carmichael Numbers')
    plt.ylabel('Rate of Fooling the Test')
    plt.title('Success Rate of Carmichael Numbers in Fooling Miller-Rabin Test')
    plt.xticks(x + width*1.5, carmichael_numbers)
    plt.legend()
    plt.savefig('carmichael_success_rate.png')
    plt.close()

def print_summary(all_results: List[Dict], carmichael_results: Dict[int, Dict[int, float]]) -> None:
    """
    Print a comprehensive summary of the experiments.
    """
    print("\nSUMMARY OF MILLER-RABIN TEST EXPERIMENTS")
    print("=" * 60)
    
    # Overall results for each k
    print("\nOverall Results by Number of Rounds:")
    for result in all_results:
        k = result['k']
        print(f"\nk = {k}:")
        print(f"  Total composites tested: {result['total_composites']:,}")
        print(f"  False positives found: {result['total_false_positives']:,}")
        print(f"  Error rate: {result['overall_error_rate']:.8f}")
    
    # Carmichael numbers analysis
    print("\nCarmichael Numbers Analysis:")
    carmichael_numbers = sorted(list(carmichael_results[list(carmichael_results.keys())[0]].keys()))
    
    for n in carmichael_numbers:
        print(f"\nNumber {n}:")
        for k in sorted(carmichael_results.keys()):
            success_rate = carmichael_results[k][n]
            print(f"  k={k}: fooled the test {success_rate*100:.1f}% of the time")
    
    # Best false positives
    print("\nMost Successful False Positives:")
    for result in all_results:
        k = result['k']
        print(f"\nk = {k}:")
        all_fps = []
        for range_data in result['range_data']:
            all_fps.extend(range_data['false_positives'])
        if all_fps:
            print("  Top 5 false positives:", ", ".join(str(n) for n in sorted(all_fps)[:5]))
        else:
            print("  No false positives found")

if __name__ == "__main__":
    # Generate complex ranges for testing
    ranges = generate_complex_ranges()
    
    # Test with various k values
    k_values = [1, 2, 3, 4, 5, 10] #15, 20, 25, 30, 40, 50]
    
    # Run experiments for different k values
    print("Starting Miller-Rabin experiments...")
    all_results = []
    carmichael_results = {}
    
    for k in k_values:
        # Run main experiment
        result = run_experiment(ranges, k)
        all_results.append(result)
        
        # Test Carmichael numbers
        carmichael_results[k] = test_carmichael_numbers(k)
    
    # Create visualizations
    print("\nGenerating visualization plots...")
    create_visualizations(all_results, carmichael_results)
    
    # Save detailed results to file
    save_results(all_results, carmichael_results)
    
    # Print brief summary to console
    print("\nExperiment completed!")
    print(f"Tested {len(k_values)} different k values: {k_values}")
    print(f"Total ranges tested: {len(ranges)}")
    print(f"Number of runs per range: 1000")
    print("See the output file for detailed results.") 