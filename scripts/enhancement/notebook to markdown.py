import nbformat
from nbconvert import MarkdownExporter
import os


def notebook_to_markdown(notebook_path, output_dir=None):
    # Load the notebook file
    with open(notebook_path, 'r', encoding='utf-8') as file:
        notebook_content = nbformat.read(file, as_version=4)

    # Initialize the Markdown exporter
    markdown_exporter = MarkdownExporter()

    # Convert the notebook to Markdown
    (body, resources) = markdown_exporter.from_notebook_node(notebook_content)

    # Define the output path
    if output_dir is None:
        output_dir = os.path.dirname(notebook_path)
    markdown_path = os.path.join(output_dir, os.path.splitext(os.path.basename(notebook_path))[0] + '.md')

    # Write the Markdown content to a file
    with open(markdown_path, 'w', encoding='utf-8') as file:
        file.write(body)

    print(f"Markdown file created at: {markdown_path}")
import os

class FileSystemAccess:
    def __init__(self, base_path=None):
        if base_path and os.path.exists(base_path):
            os.chdir(base_path)
        self.base_path = os.getcwd()

    def go_to_directory(self, path):
        if os.path.exists(path):
            os.chdir(path)
            self.base_path = os.getcwd()
        else:
            raise FileNotFoundError(f"Directory '{path}' does not exist.")

    def make_directory(self, path):
        os.makedirs(path, exist_ok=True)

    def get_all_files(self, directory=None):
        if directory is None:
            directory = self.base_path
        return [os.path.join(directory, file) for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]

    def get_files_with_extension(self, extension, directory=None):
        if directory is None:
            directory = self.base_path
        return [os.path.join(directory, file) for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file)) and file.endswith(extension)]

# Example usage
fs = FileSystemAccess()

# Go to a specific directory
fs.go_to_directory('/home/kid-a/Documents/projects/num-rel-hub/notebooks/educative/nrpy-tutorial/')


# Get all .txt files in the current directory
notebooks = fs.get_files_with_extension('.ipynb')

print(notebooks)

for notebook in notebooks:
    output_dir = "/home/kid-a/Documents/projects/num-rel-hub/docs/notes/"  # Optional, defaults to the same directory as the notebook
    notebook_to_markdown(notebook, output_dir)
