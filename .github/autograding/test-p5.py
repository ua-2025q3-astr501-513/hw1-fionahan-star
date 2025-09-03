# This is a pytest script to test the YAML file submission for
# assignment 5.
# It ensures that the file is valid, properly formatted, and meets the
# requirements.

import yaml

def test_format():
    """
    Test the format of the YAML file 'phys305_hw1/a5.yaml' to ensure:
    1. It is a valid YAML file.
    2. It contains exactly three project ideas.
    3. Each project idea is a dictionary with the required keys: 'topic' and 'idea'.
    """

    # Open and read the YAML file
    with open("hw1/p5.yaml") as f:
        try:
            # Parse the YAML content safely
            l = yaml.safe_load(f)
        except yaml.YAMLError as e:
            # Print the error if YAML parsing fails
            print("YAML Parsing Error:", e)
            assert False, "Invalid YAML file format"

        # Ensure the loaded content is a list and contains exactly three project ideas
        assert isinstance(l, list), "The YAML file must contain a list of project ideas."
        assert len(l) == 3, f"Expected exactly 3 project ideas, but found {len(l)}."

        # Iterate through the list and check that each dictionary has the required keys
        for d in l:
            assert isinstance(d, dict), "Each project entry must be a dictionary."
            assert all(k in d for k in ['topic', 'idea']), (
                f"Each project must contain 'topic' and 'idea' keys. Missing in entry: {d}"
            )

            # Print the validated dictionary entry
            print("Validated project entry:", d)
