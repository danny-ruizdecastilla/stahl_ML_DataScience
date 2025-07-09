import os
import re

def update_com_file(file_path, new_theory):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    updated_lines = []
    for line in lines:
        if line.strip().startswith("#"):
            # Replace the existing functional/basis at the beginning of the line
            line = re.sub(r'^#\s*\S+/\S+', f"# {new_theory}", line)
        updated_lines.append(line)

    with open(file_path, 'w') as f:
        f.writelines(updated_lines)

    print(f"Updated: {file_path}")

def main():
    directory = os.getcwd()
    print(f"Searching for .com files in: {directory}")
    new_theory = input("Enter new level of theory (e.g., B3LYP/6-311+G(2d,p)): ").strip()
    if not re.match(r'^\w+/\S+$', new_theory):
        print("❌ Invalid format. Use something like: B3LYP/6-311+G(2d,p)")
        return

    com_files = [f for f in os.listdir(directory) if f.endswith(".com")]
    if not com_files:
        print("No .com files found in the directory.")
        return

    for com_file in com_files:
        update_com_file(os.path.join(directory, com_file), new_theory)

if __name__ == "__main__":
    main()

