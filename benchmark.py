import random
import sys

def generate_guide(guide_length, pam_len, mismatches):
    letters = ['A','C','G','T']
    guide = "".join([random.choice(letters) for i in range(guide_length)])
    pam = "N"*pam_len
    return f"{guide}{pam} {mismatches}"

def generate_input(num_guides, guide_length, pam, mismatches, output_file):
    endl = "\n"
    guides = [generate_guide(guide_length, len(pam), mismatches) for i in range(num_guides)]
    content = f"""{'N'*guide_length}{pam}
{endl.join(guides)}
"""
    # Save to file
    with open(output_file, 'w') as f:
        f.write(content)
    return content

# Example usage:
if __name__ == "__main__":
    output_file = "input.txt"
    generate_input(20, 20, "NGG", 2, output_file)
    """num_guides = 20
    guide_length = 20
    pam = "NGG"
    mismatches = 2
    generate_input(num_guides, guide_length, pam, mismatches, output_file)"""