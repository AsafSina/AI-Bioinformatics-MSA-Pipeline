import os
import subprocess
import sys
from math import log2
from collections import Counter

# --- MSA EXECUTION FUNCTIONS ---

def run_tcoffee(input_path, output_path):
    """Executes the T-Coffee Multiple Sequence Aligner."""
    print(f"\n🔧 Running T-Coffee on: {os.path.basename(input_path)}")
    try:
        subprocess.run([
            "t_coffee", "-seq", input_path,
            "-method=t_coffee_msa",
            "-output", "fasta_aln",
            "-outfile", output_path,
            "-multi_core", "8"
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"✅ T-Coffee output saved to: {output_path}")
        return output_path
    except subprocess.CalledProcessError as e:
        print(f"❌ T-Coffee alignment failed: {e.stderr.decode().strip()}")
    except Exception as e:
        print(f"❌ T-Coffee unexpected error: {e}")
    return None

def run_mafft(input_path, output_path):
    """Executes the MAFFT Multiple Sequence Aligner."""
    print(f"\n🔧 Running MAFFT on: {os.path.basename(input_path)}")
    try:
        with open(output_path, "w") as out_f:
            subprocess.run([
                "mafft", "--auto", input_path
            ], stdout=out_f, stderr=subprocess.PIPE, check=True)
        print(f"✅ MAFFT output saved to: {output_path}")
        return output_path
    except subprocess.CalledProcessError as e:
        print(f"❌ MAFFT alignment failed: {e.stderr.decode().strip()}")
    except Exception as e:
        print(f"❌ MAFFT unexpected error: {e}")
    return None

def run_clustalo(input_path, output_path):
    """Executes the Clustal Omega Multiple Sequence Aligner."""
    print(f"\n🔧 Running Clustal Omega on: {os.path.basename(input_path)}")
    try:
        subprocess.run([
            "clustalo", "-i", input_path, "-o", output_path,
            "--force", "--threads=8"
        ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"✅ Clustal Omega output saved to: {output_path}")
        return output_path
    except subprocess.CalledProcessError as e:
        print(f"❌ Clustal Omega alignment failed: {e.stderr.decode().strip()}")
    except Exception as e:
        print(f"❌ Clustal Omega unexpected error: {e}")
    return None

# --- SCORING FUNCTIONS ---

def read_fasta(file_path):
    """Reads a FASTA MSA file and returns a list of sequences."""
    sequences = []
    seq = ""
    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if seq:
                        sequences.append(seq)
                        seq = ""
                else:
                    seq += line.replace(' ', '')
            if seq:
                sequences.append(seq)
        return sequences
    except Exception as e:
        print(f"❌ FASTA read error: {e}")
        return []

def shannon_entropy(column):
    """Calculates Shannon entropy for one MSA column."""
    counts = Counter(column)
    total = len(column)
    entropy = 0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * log2(p)
    return entropy

def calculate_scores(sequences):
    """Computes gap %, match % and entropy for the whole MSA."""
    if not sequences:
        return 0, 0, 0

    seq_len = len(sequences[0])
    num_seqs = len(sequences)

    gap_count = 0
    match_count = 0
    entropy_sum = 0

    for col_idx in range(seq_len):
        column = [seq[col_idx] for seq in sequences]
        entropy_sum += shannon_entropy(column)

        # Gaps
        gap_count += column.count('-')

        # Matches (most common AA)
        most_common = max(set(column), key=column.count)
        match_count += column.count(most_common)

    total_positions = seq_len * num_seqs
    gap_percentage = (gap_count / total_positions) * 100 if total_positions else 0
    match_percentage = (match_count / total_positions) * 100 if total_positions else 0
    avg_entropy = entropy_sum / seq_len if seq_len else 0

    return gap_percentage, match_percentage, avg_entropy

def score_and_report(file_path):
    """Scores the alignment file and prints the results."""
    sequences = read_fasta(file_path)
    if sequences:
        gap_p, match_p, entropy = calculate_scores(sequences)
        print(f"\n--- 📊 Scoring Results ({os.path.basename(file_path)}) ---")
        print(f"  Gap Percentage: {gap_p:.3f}%")
        print(f"  Match Percentage: {match_p:.3f}%")
        print(f"  Average Shannon Entropy: {entropy:.4f}")
        return {
            "file": os.path.basename(file_path),
            "gap_p": gap_p,
            "match_p": match_p,
            "entropy": entropy
        }
    return None

# --- MAIN MENU PIPELINE (Integrated) ---

def run_msa_pipeline_menu():
    print("--- 🧬 FULL MSA Pipeline (Alignment and Scoring) Initiated ---")

    input_path = input("📂 Enter the path to a FASTA file or directory: ").strip()
    if not os.path.exists(input_path):
        print(f"❌ Path not found: {input_path}")
        sys.exit(1)

    output_dir = input("📁 Enter the directory to save the outputs: ").strip()
    os.makedirs(output_dir, exist_ok=True)

    print("\n--- Alignment Tool Selection ---")
    print("1: MAFFT")
    print("2: Clustal Omega")
    print("3: T-Coffee")
    print("4: Run All Aligners (For Comparative Analysis)")
    choice = input("👉 Select the MSA tool(s) to run (1-4): ").strip()

    # Determine input files
    files = []
    if os.path.isfile(input_path):
        files = [input_path]
    elif os.path.isdir(input_path):
        files = [os.path.join(input_path, f) for f in os.listdir(input_path) if f.endswith((".fasta", ".fa", ".fna"))]
        if not files:
            print("❌ No FASTA files found in the directory.")
            sys.exit(1)

    # Run alignment and scoring for each file
    global_report = []
    for f in files:
        basename = os.path.splitext(os.path.basename(f))[0]
        output_paths = []

        def run_tool(tool_func, tool_name):
            out_path = os.path.join(output_dir, f"{basename}_{tool_name}.fasta")
            result_path = tool_func(f, out_path)
            if result_path:
                output_paths.append(result_path)

        if choice == '1':
            run_tool(run_mafft, 'mafft')
        elif choice == '2':
            run_tool(run_clustalo, 'clustal')
        elif choice == '3':
            run_tool(run_tcoffee, 'tcoffee')
        elif choice == '4':
            run_tool(run_tcoffee, 'tcoffee')
            run_tool(run_mafft, 'mafft')
            run_tool(run_clustalo, 'clustal')
        else:
            print("🛑 Invalid selection. Pipeline aborted.")
            sys.exit(1)

        # Score successful alignments
        for out_f in output_paths:
            score_data = score_and_report(out_f)
            if score_data:
                global_report.append(score_data)

    # Print the Global Comparative Report (Meaningful for choice 4)
    if choice == '4' and global_report:
        print("\n" + "="*70)
        print("🏆 COMPARATIVE ANALYSIS SUMMARY")
        print("="*70)
        print(f"{'Alignment File':<30} | {'Gap (%)':>10} | {'Match (%)':>10} | {'Entropy':>10} | {'Tool':<10}")
        print("-" * 70)
        for data in global_report:
            # Extract tool name from filename for better report
            tool_name = data['file'].split('_')[-1].split('.')[0]
            print(f"{data['file']:<30} | {data['gap_p']:>10.3f} | {data['match_p']:>10.3f} | {data['entropy']:>10.4f} | {tool_name.capitalize():<10}")
        print("-" * 70)
        print("💡 Note: Lower Entropy indicates higher alignment conservation.")


if __name__ == "__main__":
    run_msa_pipeline_menu()