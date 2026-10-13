
import subprocess
import re

filename = "testC.in"
N = 50

# Read the original input file
with open(filename, "r") as f:
    lines = f.readlines()

start = lines.index("#card 2: dx\n") + 1
end = lines.index("#card 3: number of Discrete ordinates (DO) \n")

def calculate_keff(a):
    dx = a / N

    new_lines = lines[:start]
    for i in range(0, N, 5):
        new_lines.append((f"{dx:.12f} " * 5) + "\n")
    new_lines += lines[end:]

    with open(filename, "w") as f:
        f.writelines(new_lines)

    # Run solver
    result = subprocess.run(
        ["../eigensolver/snif", filename],
        capture_output=True,
        text=True,
        check=True
    )

    # Extract Keff
    keff = float(result.stdout.splitlines()[12].split()[-1])

    return keff


low = 5.0
high = 10.0

for i in range(50):
    a = (low + high) / 2
    keff = calculate_keff(a)

    print(f"{i+1}: a={a:.12f}, keff={keff:.12f}")

    if abs(keff - 1.0) < 1e-12:
        break

    if keff < 1.0:
        low = a
    else:
        high = a

print(f"\nCritical width = {a:.12f} cm")
print(f"Keff = {keff:.12f}")

