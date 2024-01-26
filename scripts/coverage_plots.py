try:
    import matplotlib.pyplot as plt
except:
    print("pyplot import failed")
try:
    import os
except:
    print("OS import failed")
try:
    import numpy as np
except:
    print("numpy import failed")
try:
    from sklearn.preprocessing import PolynomialFeatures
except:
    print("sklearn.preprocessing import failed")
try:
    from sklearn.linear_model import LinearRegression
except:
    print("sklearn.linear_model import failed")

directory_path = snakemake.output['path']

if not os.path.exists(directory_path):
    try:
        os.mkdir(directory_path)
        print(f"Directory '{directory_path}' created successfully.")
    except OSError as error:
        print(f"Failed to create directory '{directory_path}': {error}")
else:
    print(f"Directory '{directory_path}' already exists.")
          
for scaffold in snakemake.input:
    print("Now processing " + scaffold)
    # Read data from the file
    base_positions = []
    frequencies = []
    
    with open(scaffold, "r") as f:
        for line in f:
            parts = line.strip().split()
            base_positions.append(int(parts[1]))  # Assuming base positions are in the third column
            frequencies.append(int(parts[2]))     # Assuming frequencies are in the second column
    
    # Convert lists to numpy arrays
    base_positions_np = np.array(base_positions).reshape(-1, 1)  # Reshape to column vector
    frequencies_np = np.array(frequencies)
    
    # Fit polynomial regression model
    poly_features = PolynomialFeatures(degree=10)  # Adjust the degree as needed
    base_positions_poly = poly_features.fit_transform(base_positions_np)
    model = LinearRegression()
    model.fit(base_positions_poly, frequencies_np)
    
    # Predict frequencies using the model
    base_positions_range = np.linspace(min(base_positions), max(base_positions), 100).reshape(-1, 1)
    base_positions_range_poly = poly_features.transform(base_positions_range)
    predicted_frequencies = model.predict(base_positions_range_poly)

    # Compute moving average
    window_size = 10000  # Adjust the window size as needed
    moving_avg = np.convolve(frequencies, np.ones(window_size) / window_size, mode='valid')
    
    # Adjust base positions to match the moving average length
    adjusted_base_positions = base_positions[window_size-1:]
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    

    # Create the plot
    plt.figure(figsize=(13, 6))  # Adjust the figure size as needed
    plt.plot(base_positions, frequencies, color='lightblue', linewidth=1.25, label='Coverage depth')
    plt.plot(base_positions_range, predicted_frequencies, color='blue', linestyle='--', label='Trend Line')
    plt.plot(adjusted_base_positions, moving_avg, color='red', linewidth=.75, label='Moving Average (10,000 bp)')

    # Add labels and title
    plt.xlabel('Base position along sequence')
    plt.ylabel('Depth of coverage')
    plt.title('Depth of coverage vs. Base position')
   
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.7)
   
    # Show legend
    plt.legend()

    # Save the plot to file
    out_png = snakemake.output['path'] + "/" + scaffold.split("/")[-1] + ".png"
    plt.savefig(out_png)  # Save as PNG file