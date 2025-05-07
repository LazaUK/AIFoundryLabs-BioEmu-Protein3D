# Gradio Web UI app to predict 3D protein structures with BioEmu (Biomolecular Emulator)
This demo app provides a user-friendly Web UI wrapper for predicting 3D protein structures with the [BioEmu](https://github.com/microsoft/bioemu) AI model, released as open source by the Microsoft Research team on [Azure AI Foundry Labs](https://ai.azure.com/labs).

With this app, you can:
- Input protein sequences in plain text or FASTA format,
- Visualise generated protein structures in an interactive 3D viewer,
- Download generated protein structures.

## Table of contents:
- [Part 1: Local Use of Companion Docker Image]()
- [Part 2: Cloud Deployment to Azure Web App]()
- [Part 3: User Experience - Gradio UI]()
- [Demo video on YouTube]()
- [Acknowledgments and Citation]()

## Part 1: Local Use of Companion Docker Image
> [!NOTE]
> This repo comes with a pre-built Docker image, available at: `ghcr.io/lazauk/bioemu-webapp:latest`.

1. To pull the image locally, use the Docker command below. Alternatively, if you are deploying to an Azure Web app, skip this and move to Part 2.
``` PowerShell
docker pull ghcr.io/lazauk/bioemu-webapp:latest
```
2. You can start the Docker container with the following command:
``` PowerSehll
docker run -p 7860:7860 ghcr.io/lazauk/bioemu-webapp:latest
```
3. The Web app can now be accessed at http://127.0.0.1:7860.

## Part 2: Cloud Deployment to Azure Web App
See our [detailed deployment guide](DEPLOYMENT.md) for instructions on deploying this application to Azure Web App for Containers.

> [!WARNING]
> - First-time protein structure generation may be slow as the BioEmu model is downloaded from HuggingFace
> - Larger proteins require more computation time
> - For optimal performance, we recommend deploying with at least 4GB of RAM

## Part 3: User Experience - Gradio UI
1. Enter a protein sequence or choose from sample sequences
2. Set the number of structures to generate and other parameters
3. Click "Generate Structures"
4. Explore the 3D visualizations and analysis
5. Download the resulting structures if desired

## Demo video on YouTube
Visit our demo application at [https://bioemu-app-service.azurewebsites.net](https://bioemu-app-service.azurewebsites.net)

## Acknowledgments and Citation
- [Microsoft BioEmu](https://github.com/microsoft/bioemu) - The core AI model
- [Gradio](https://gradio.app/) - The web UI framework

If you use this demo or the AI model itself, please cite the original BioEmu paper of Microsoft Research team:
``` bibtex
@article {BioEmu2024,
    author = {Lewis, Sarah and Hempel, Tim and Jim{\'e}nez-Luna, Jos{\'e} and Gastegger, Michael and Xie, Yu and Foong, Andrew Y. K. and Satorras, Victor Garc{\'\i}a and Abdin, Osama and Veeling, Bastiaan S. and Zaporozhets, Iryna and Chen, Yaoyi and Yang, Soojung and Schneuing, Arne and Nigam, Jigyasa and Barbero, Federico and Stimper, Vincent and Campbell, Andrew and Yim, Jason and Lienen, Marten and Shi, Yu and Zheng, Shuxin and Schulz, Hannes and Munir, Usman and Clementi, Cecilia and No{\'e}, Frank},
    title = {Scalable emulation of protein equilibrium ensembles with generative deep learning},
    year = {2024},
    doi = {10.1101/2024.12.05.626885},
    journal = {bioRxiv}
}
```
