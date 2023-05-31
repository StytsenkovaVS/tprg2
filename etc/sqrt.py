import numpy as np

def triangular_transform(data):
    """
    Transform a pseudorandom number sequence into another pseudorandom number
    sequence with triangular distribution.
    
    Parameters:
    data (numpy.ndarray): An array of pseudorandom numbers with a uniform distribution
    
    Returns:
    numpy.ndarray: An array of pseudorandom numbers with a triangular distribution
    """
    transformed = np.sqrt(data)
    return transformed

# Generate 1000 random numbers from a uniform distribution between 0 and 1
data = np.random.uniform(size=10000)

# Transform the data using the triangular_transform function
transformed = triangular_transform(data)

# Plot a histogram of the transformed data
import matplotlib.pyplot as plt
plt.hist(transformed, bins=20)

# Add labels and a title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Triangular Distribution')

# Show the plot
plt.show()
