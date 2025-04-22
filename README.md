# Azure AI Foundry Labs: Gradio Web UI app to demo Biomolecular Emulator (BioEmu)
This application provides a user-friendly Web interface for interacting with the [BioEmu](https://github.com/microsoft/bioemu) protein structure prediction model. With this demo, you can:
- Input protein sequences in plain text or FASTA format
- Visualise generated protein structures in an interactive 3D viewer
- Download generated protein structures



## Using Companion Docker image
This repository is set up to publish Docker images to GitHub Container Registry. The pre-built image is available at:

```
ghcr.io/GITHUB_USERNAME/bioemu-demo:latest
```

To pull the image:

```bash
docker pull ghcr.io/GITHUB_USERNAME/bioemu-demo:latest
```

## Deploying to Azure Web App
See our [detailed deployment guide](DEPLOYMENT.md) for instructions on deploying this application to Azure Web App for Containers.

> NOTE
> - First-time protein structure generation may be slow as the BioEmu model is downloaded from HuggingFace
> - Larger proteins require more computation time
> - For optimal performance, we recommend deploying with at least 4GB of RAM

## Usage
1. Enter a protein sequence or choose from sample sequences
2. Set the number of structures to generate and other parameters
3. Click "Generate Structures"
4. Explore the 3D visualizations and analysis
5. Download the resulting structures if desired

## Live Demo
Visit our demo application at [https://bioemu-app-service.azurewebsites.net](https://bioemu-app-service.azurewebsites.net)

## Acknowledgments
- [Microsoft BioEmu](https://github.com/microsoft/bioemu) - The core AI model
- [Gradio](https://gradio.app/) - The web UI framework

## Citation
If you use this demo in your research, please cite the original BioEmu paper of Microsoft Research team:

```bibtex
@article {BioEmu2024,
    author = {Lewis, Sarah and Hempel, Tim and Jim{\'e}nez-Luna, Jos{\'e} and Gastegger, Michael and Xie, Yu and Foong, Andrew Y. K. and Satorras, Victor Garc{\'\i}a and Abdin, Osama and Veeling, Bastiaan S. and Zaporozhets, Iryna and Chen, Yaoyi and Yang, Soojung and Schneuing, Arne and Nigam, Jigyasa and Barbero, Federico and Stimper, Vincent and Campbell, Andrew and Yim, Jason and Lienen, Marten and Shi, Yu and Zheng, Shuxin and Schulz, Hannes and Munir, Usman and Clementi, Cecilia and No{\'e}, Frank},
    title = {Scalable emulation of protein equilibrium ensembles with generative deep learning},
    year = {2024},
    doi = {10.1101/2024.12.05.626885},
    journal = {bioRxiv}
}
```
