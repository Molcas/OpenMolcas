import random
from itertools import combinations

def create_site_dictionary(initial_key):
    """
    Create all possible dictionaries (35) based on different ways to split sites into groups.
    
    Args:
        initial_key: An integer key that determines the site that must be in group 1
        
    Returns:
        List of dictionaries with keys 1-40 and values assigned according to rules
    """
    # Create empty list to store all dictionaries
    all_dicts = []
    
    # Map keys to their respective sites
    sites = {}
    for site_num in range(1, 9):
        start_key = (site_num - 1) * 5 + 1
        sites[site_num] = list(range(start_key, start_key + 5))
    
    # Determine which site the initial_key belongs to
    initial_key_site = None
    for site_num, site_keys in sites.items():
        if initial_key in site_keys:
            initial_key_site = site_num
            break
    
    if initial_key_site is None:
        raise ValueError(f"Initial key {initial_key} is not in range 1-40")
    
    # Get all other sites
    other_sites = [site for site in range(1, 9) if site != initial_key_site]
    
    # Generate all combinations of 3 sites from the remaining 7 sites
    # Each combination plus the initial_key_site will form group 1
    for combo in combinations(other_sites, 3):
        # Create new dictionary
        result_dict = {}
        
        # Group 1 consists of initial_key_site and 3 other sites
        group1_sites = [initial_key_site] + list(combo)
        group2_sites = [site for site in range(1, 9) if site not in group1_sites]
        
        # Set the first key of each site to value 3
        for site_keys in sites.values():
            result_dict[site_keys[0]] = 3
        
        # Since initial_key is in group1_sites by definition
        group1_value, group2_value = 1, 2
        
        # Assign remaining values based on group membership
        for site_num in group1_sites:
            for key in sites[site_num][1:]:  # Skip first key in each site
                result_dict[key] = group1_value
                
        for site_num in group2_sites:
            for key in sites[site_num][1:]:  # Skip first key in each site
                result_dict[key] = group2_value
        
        all_dicts.append(result_dict)
    
    return all_dicts

def main():
    # Example usage
    initial_key = 1  # This can be any key from 1 to 40
    site_dicts = create_site_dictionary(initial_key)
    
    # Print information about the generated dictionaries
    print(f"Generated {len(site_dicts)} different dictionaries")
    print(f"First dictionary example:")
    for key in sorted(site_dicts[1].keys()):
        print(f"Key: {key}, Value: {site_dicts[1][key]}")

if __name__ == "__main__":
    main()
