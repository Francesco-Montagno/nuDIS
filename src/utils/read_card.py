import os
import sys

class InputCard:
    def __init__(self, card_path):
        """
        Initializes the InputCard reader.
        
        :param card_path: Path to the configuration file (e.g., 'card/run_card.dat')
        """
        self.card_path = card_path
        self.params = {}
        
        if not os.path.exists(card_path):
            print(f"CRITICAL ERROR: The card file '{card_path}' was not found.")
            sys.exit(1)
            
        self._parse_card()

    def _parse_card(self):
        """Reads the file line by line, ignoring comments and whitespace."""
        with open(self.card_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                raw_line = line.strip()
                
                # Skip empty lines or full-line comments
                if not raw_line or raw_line.startswith('#'):
                    continue
                
                # Remove inline comments (e.g., var = 1 # comment)
                line_content = raw_line.split('#')[0].strip()
                
                if '=' in line_content:
                    key, val = line_content.split('=', 1)
                    key = key.strip()
                    val = val.strip()
                    self.params[key] = self._convert_type(val)
                else:
                    print(f"Warning: Line {line_num} ignored (incorrect format): {raw_line}")

    def _convert_type(self, value):
        """Attempts to convert the string value to int or float."""
        # Try converting to integer
        try:
            return int(value)
        except ValueError:
            pass
            
        # Try converting to float
        try:
            return float(value)
        except ValueError:
            pass
            
        # Fallback: Return as string
        return value

    def get(self, key):
        """
        Retrieves a value from the parameters. 
        Crashes the program if the key is missing to prevent silent errors.
        """
        if key not in self.params:
            print(f"ERROR: Parameter '{key}' is missing in the card '{self.card_path}'")
            sys.exit(1)
        return self.params[key]