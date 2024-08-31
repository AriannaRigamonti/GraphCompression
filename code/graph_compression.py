#!/usr/bin/env python3

import argparse
import requests
import os 
import pandas as pd
import io 
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering

# Define global paths to save all the following outputs
OUTPUT_FOLDER = 'output_files'

def create_directory(path): # check if the folder exists, otherwise create it. This code will be used many times so I decided to create a function. 
    if not os.path.exists(path):
        os.makedirs(path)

def get_unique_filepath(filepath):
    """Generate a unique filepath by appending a suffix if the file already exists."""
    base, extension = os.path.splitext(filepath) # take file name
    counter = 1
    while os.path.exists(filepath):
        filepath = f"{base}_{counter}{extension}"
        counter += 1
    return filepath

class Protein:
    def __init__(self, overwrite, name, species, nodes_number):
        self.overwrite = overwrite # default = False
        self.name = name
        self.species = species # default = 9606 (human)
        self.nodes_number = nodes_number
        self.identifier = self.get_identifier() # a method of this class
        self.interaction_network = self.get_interaction_network() # a method of this class
        
    def get_identifier(self): 
        """Maps common protein names, synonyms and UniProt identifiers into STRING identifiers, that will be used to retreive the interaction network in a unique way.
        The method sends a GET request to the STRING API to retrieve the protein identifier and returns it if the request is successful. 
        If the request fails, it prints an error message and returns None."""

        if self.name: # if protein name was given we will use this method (if instead dataset was given, we skip this part)
            print("Getting protein identifier...")

            # URL of the network in the STRING API
            string_identifier_url = "https://string-db.org/api/tsv/get_string_ids?"
            
            # Parameters for the request
            parameter_identifier = {
                "identifier": self.name,
                "species": self.species #NCBI taxon identifiers (e.g. Human is 9606, see: STRING organisms).
                }
            
            # Send GET request to STRING API
            http_response = requests.get(string_identifier_url, params=parameter_identifier)
            
            # Check the status of the request
            if http_response.status_code == 200:

                # Read the response into a pd DataFrame
                df_identifier = pd.read_csv(io.StringIO(http_response.text), sep='\t')
                
                # Get the value from the second column
                string_identifier = df_identifier.iloc[0, 1]

            else:
                # If the request fails, print an error message
                print(f"Error in request: {http_response.status_code}")
            
            return string_identifier # it will be used to get the interactoin network when only protein name is given
        
    def get_interaction_network(self): 
        # Retrieves the network interactions for the input protein 
        print("Getting interaction network...")
    
        # URL of the interaction partners in the STRING API
        string_interaction_network = "https://string-db.org/api/tsv/network?"

        # Parameters for the request
        parameters_partners = {
            "identifiers": self.identifier,
            "add_nodes": self.nodes_number,

        }

        # Send GET request to STRING API
        http_response = requests.get(string_interaction_network, params=parameters_partners)

        # Check the status of the request
        if http_response.status_code == 200:

            # Parse the response as a DataFrame
            interaction_network = pd.read_csv(io.StringIO(http_response.text), sep='\t')

            # Save the file 
            subdirectory = os.path.join(OUTPUT_FOLDER, 'interaction_network', self.name)
            create_directory(subdirectory)

            csv_file_path = os.path.join(subdirectory, f"{self.name}_{self.nodes_number}_nodes_interaction_network.csv")
                           
            if self.overwrite is False:   
                csv_file_path = get_unique_filepath(csv_file_path) # do not overwrite, but add _counter at the end 
            
            interaction_network.to_csv(csv_file_path, index=False)

            return interaction_network

        else:
            print(f"Error > Status code: {http_response.status_code}")

            return None
    
class DataFrame:
    def __init__(self, input_df):
        self.input_df = pd.read_csv(input_df)
        self.df = input_df
        self.input_df_name =  self.get_input_df_name() # a method of this class

    def get_input_df_name(self):
        """Returns the name of the input DataFrame file."""
        input_df_name = os.path.basename(self.df)
        return input_df_name 

class Network:
    def __init__(self, overwrite, cluster_number, protein = None, dataframe = None):
        self.overwrite = overwrite
        self.cluster_number = cluster_number
        if protein is not None: 
            self.protein = protein
            self.dataframe = None
        elif dataframe is not None: 
            self.dataframe = dataframe
            self.protein = None
        self.G = self.create_graph() # a method of this class
        self.compressed_graph, self.node_original_mapping = self.compress_graph() # a method of this class
        self.compressed_interaction_network = self.get_compressed_interaction_network() # a method of this class
        
    def create_graph(self):
        """Creates a graph from the interaction network data"""
        print("Creating graph...") 

        # if dataframe is given it creates the network directly (otherwise the interaction data is retrieved from the protein name).
        if self.protein is not None:
            interaction_network = self.protein.interaction_network
        else:
            interaction_network = self.dataframe.input_df

        # Create an empty graph
        G = nx.Graph()

        try:
            # Add nodes and edges to the graph
            for index, row in interaction_network.iterrows():
                G.add_node(row['preferredName_A'])
                G.add_node(row['preferredName_B'])
                G.add_edge(row['preferredName_A'], row['preferredName_B'], weight=row['score'])

                # note: score is the combined score of: 
                        # nscore: gene neighborhood score
                        # fscore: gene fusion score
                        # pscore: phylogenetic profile score
                        # ascore: coexpression score
                        # escore: experimental score
                        # dscore: database score
                        # tscore: textmining score

            # Draw the graph
            plt.figure(figsize=(12, 8))
            pos = nx.spring_layout(G, seed=42)  # Nodes position
            nx.draw(G, pos, with_labels=True, node_size=1000, node_color='skyblue', font_size=8, font_weight='bold') 
            plt.gcf().suptitle(f'Protein Interaction Network ({len(G.nodes)-1} nodes) for {self.protein.name}') if self.protein else plt.gcf().suptitle(f'Protein Interaction Network ({len(G.nodes)-1} nodes) for {self.dataframe.input_df_name}')

            # Store the graph png in the directory
            subdirectory = os.path.join(OUTPUT_FOLDER, 'original_graph', self.protein.name if self.protein else self.dataframe.input_df_name)
            create_directory(subdirectory)
            if self.protein:
                file_path = os.path.join(subdirectory, f'{self.protein.name if self.protein else self.dataframe.input_df_name}_{len(G.nodes)-1}_nodes_original_graph.png')
            else:
                file_path = os.path.join(subdirectory, f'{self.dataframe.input_df_name}_original_graph.png')

            if self.overwrite is not False:
                file_path = get_unique_filepath(file_path) # do not overwrite
            
            plt.savefig(file_path, format='png')

        except (pd.errors.EmptyDataError, FileNotFoundError) as e:
            print(f"Error creating graph: {e}")
            
        return G
        
    def cluster_nodes(self, num_clusters):
        """Clusters the nodes of the graph using Spectral Clustering."""

        try:
            # Convert the graph to a numpy adjacency matrix
            adjacency_matrix = nx.to_numpy_array(self.G)

            # Perform Spectral Clustering (this is one of the many options to compress graphs. I choose this just because I read that it is widely used with networks).
            spectral = SpectralClustering(n_clusters=num_clusters, affinity='precomputed', random_state=42)
            cluster_labels = spectral.fit_predict(adjacency_matrix)

            # Map nodes to their respective clusters
            node_cluster_mapping = {node: cluster_labels[i] for i, node in enumerate(self.G.nodes())}
            
            return node_cluster_mapping

        except Exception as e:
            print(f"Error in clustering nodes: {e}")
            return None
    
    def contract_edges(self, node_cluster_mapping):
        """Contracts the edges between nodes belonging to the same cluster."""

        # Create a new graph for the contracted version
        compressed_graph = nx.Graph()

        try:
            # Iterate through the original graph edges and add edges between clusters
            for u, v, data in self.G.edges(data=True):
                cluster_u = node_cluster_mapping[u]
                cluster_v = node_cluster_mapping[v]
                if (cluster_u, cluster_v) in compressed_graph.edges():
                    compressed_graph[cluster_u][cluster_v]['weight'] += data['weight']
                else:
                    compressed_graph.add_edge(cluster_u, cluster_v, weight=data['weight'])
            
            return compressed_graph

        except Exception as e:
            print(f"Error in contracting edges: {e}")
            return None

    def compress_graph(self):
        """Compress the graph by clustering nodes and contracting edges."""
        print("Compressing graph...")

        num_clusters = self.cluster_number
        node_cluster_mapping = self.cluster_nodes(num_clusters)
        compressed_graph = self.contract_edges(node_cluster_mapping)

        # Associate nodes of the compressed graph with original nodes
        node_original_mapping = {f"cluster {new_node}": [', '.join([f'{old_node}' for old_node, cluster in node_cluster_mapping.items() if cluster == new_node])] for new_node in compressed_graph.nodes()}

        # Reshape the data for DataFrame creation
        data = [(compressed_node, original_nodes[0]) for compressed_node, original_nodes in node_original_mapping.items()]

        # Create the DataFrame
        mapping_df = pd.DataFrame(data, columns=['New_node', 'Original_nodes'])

        # Save the mapping to a CSV file
        subdirectory = os.path.join(OUTPUT_FOLDER, 'node_mapping', self.protein.name if self.protein else self.dataframe.input_df_name)
        create_directory(subdirectory)
        if self.protein: 
            file_path = os.path.join(subdirectory, f'{self.protein.name}_{self.protein.nodes_number}_to_{self.cluster_number}_mapping.csv')
        else:
            file_path = os.path.join(subdirectory, f'{self.dataframe.input_df_name}_to_{self.cluster_number}_mapping.csv')
        
        if self.overwrite is not False:
            file_path = get_unique_filepath(file_path) # do not overwrite
        
        mapping_df.to_csv(file_path, index=False)

        # Creation of the compressed graph
        plt.figure(figsize=(12, 8))
        pos = nx.spring_layout(compressed_graph, seed=42)  
        nx.draw(compressed_graph, pos, with_labels=True, node_size=1000, node_color='skyblue', font_size=8, font_weight='bold')  
        
        if self.protein:
            plt.gcf().suptitle(f'Compressed Protein Interaction Network ({self.protein.nodes_number} to {num_clusters}) for {self.protein.name}')

        else:
            plt.gcf().suptitle(f'Compressed Protein Interaction Network (to {num_clusters}) for {self.dataframe.input_df_name}')
            
        # Add legend with mapping information from DataFrame
        for i, row in mapping_df.iterrows():
            plt.text(-1.0, -0.80 + 0.05*i, f"{row['New_node']}: {row['Original_nodes']}", fontsize=8, ha='left', va='top')
            plt.subplots_adjust(right=0.7, top=0.9)

        # Save the compressed graph as a PNG image
        subdirectory = os.path.join(OUTPUT_FOLDER, 'compressed_graph', self.protein.name if self.protein else self.dataframe.input_df_name)
        create_directory(subdirectory)
        if self.protein:
            file_path = os.path.join(subdirectory, f'{self.protein.name}_{self.protein.nodes_number}_to_{num_clusters}_compressed_graph_.png')
        else:
            file_path = os.path.join(subdirectory, f'to_{num_clusters}_compressed_graph_.png')

        if self.overwrite is not False:
            file_path = get_unique_filepath(file_path) # do not overwrite
        
        plt.savefig(file_path, format='png')
        
        return compressed_graph, node_original_mapping

    def get_compressed_interaction_network(self):
        """Retrieves the interaction network for the compressed graph."""
        print("Computing compressed network...")

        num_clusters = self.cluster_number

        interactions = []

        # Iterate through edges of the compressed graph
        for u, v, data in self.compressed_graph.edges(data=True):
            # Append cluster and score to the interactions 
            interactions.append((u, v, round(data['weight'], 3)))

        # Create DataFrame from interaction list
        compressed_interaction_df = pd.DataFrame(interactions, columns=['clusterA', 'clusterB', 'score'])
        compressed_interaction_df['preferred_nameA'] = ''
        compressed_interaction_df['preferred_nameB'] = ''

        # Iterate through the DataFrame and update the preferred cluster names
        for index, row in compressed_interaction_df.iterrows():
            clusterA = row['clusterA']
            clusterB = row['clusterB']

            preferred_nameA = self.node_original_mapping.get(f"cluster {clusterA}", "")
            preferred_nameB = self.node_original_mapping.get(f"cluster {clusterB}", "")

            compressed_interaction_df.at[index, 'preferred_nameA'] = ', '.join(preferred_nameA)
            compressed_interaction_df.at[index, 'preferred_nameB'] = ', '.join(preferred_nameB)

        new_column_order = ['clusterA', 'clusterB', 'preferred_nameA', 'preferred_nameB', 'score']
        compressed_interaction_df = compressed_interaction_df.reindex(columns=new_column_order)   

        # Save the DataFrame as a CSV file
        subdirectory = os.path.join(OUTPUT_FOLDER, 'compressed_interaction_network', self.protein.name if self.protein else self.dataframe.input_df_name)
        create_directory(subdirectory)
        if self.protein:
         file_path = os.path.join(subdirectory, f'{self.protein.name}_{self.protein.nodes_number}_to_{num_clusters}_compressed_interaction_network.csv')
        else:
            file_path = os.path.join(subdirectory, f'to_{num_clusters}_compressed_interaction_network.csv')
       
        if self.overwrite is not False:
            file_path = get_unique_filepath(file_path)
        
        compressed_interaction_df.to_csv(file_path, index=False)

        # Find the full path to the current script
        script_path = os.path.realpath(__file__)

        # Find the directory containing the script
        script_directory = os.path.dirname(script_path)

        print(f"Graph compression DONE!\nPlease check output_files in the folder: {script_directory}")

        return compressed_interaction_df
        
def parse_args():
    """Parse command-line arguments."""
    try: 
        parser = argparse.ArgumentParser(description="Graph Compression")
        parser.add_argument('-p', '--protein_name', type=str, required=False, help='Protein name or identifier')
        parser.add_argument('-s', '--species', type=int, default=9606, help='Taxonomy identifier of species (e.g Human is 9606)')
        parser.add_argument('-d', '--data_frame', type=str, required=False, help='Input dataframe in CSV format')
        parser.add_argument('-n', '--nodes_number', type=int, required=False, help='Number of interacting partners')
        parser.add_argument('-c', '--clusters_number', type=int, required=True, help='Number of clusters for spectral clustering')
        parser.add_argument('-o', '--overwrite', type = bool, default=False, help='Input if want to overwrite files with same name')
        args = parser.parse_args()

        # Check if both -p and -d are provided
        if args.protein_name and args.data_frame:
            parser.error("Both -p and -d options cannot be provided simultaneously. Please provide only one.")

        # Check if at least one of -p or -d is provided
        if not (args.protein_name or args.data_frame):
            parser.error("Either -p or -d option must be provided.")

        return args
    
    except Exception as e:
        print(f"An unexpected error occurred in argparse: {e}")
        exit(1)

def main():
    """Main function to run the script."""
    args = parse_args()
    protein_name = args.protein_name
    species = args.species
    input_df = args.data_frame
    nodes_number = args.nodes_number
    clusters_number = args.clusters_number
    overwrite = args.overwrite
    
    if protein_name:
        protein = Protein(overwrite, protein_name, species, nodes_number)
        network = Network(overwrite, clusters_number, protein=protein)
    elif input_df:
        df = DataFrame(input_df)
        network = Network(overwrite, clusters_number, dataframe=df)
    else:
        print("Either protein_name or data_frame must be provided.")

if __name__ == "__main__":
    main()