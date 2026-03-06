import json
import sys

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <json_file>")
    sys.exit(1)

filename = sys.argv[1]

try:
    with open(filename, 'r') as file:
        data = json.load(file)
except FileNotFoundError:
    print(f"Error: The file '{filename}' was not found.")
    sys.exit(1)
except json.JSONDecodeError:
    print(f"Error: Failed to decode JSON from '{filename}'.")
    sys.exit(1)

print("Top-level keys and their item counts:")

# Iterate through the top-level keys
for key, value in data.items():
    # Determine the count based on the type of the value
    if isinstance(value, list) or isinstance(value, dict):
        count = len(value)
        print(f"Key: '{key}', Type: {type(value).__name__}, Item Count: {count}")
    else:
        # For primitive types (string, number, boolean, null)
        count = 1 # Or some other appropriate value if you want to count individual items differently
        print(f"Key: '{key}', Type: {type(value).__name__}, Item Count: {count} (primitive value)")
