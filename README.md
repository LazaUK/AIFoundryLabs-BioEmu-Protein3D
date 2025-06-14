# Gradio Web UI app to predict 3D protein structures with BioEmu (Biomolecular Emulator)

This demo app provides a user-friendly Web UI wrapper for predicting 3D protein structures with the [BioEmu](https://github.com/microsoft/bioemu) AI model, released as open source by the Microsoft Research team on [Azure AI Foundry Labs](https://ai.azure.com/labs).

With this app, you can:
- Input protein sequences in plain text or FASTA format,
- Visualise generated protein structures in an interactive 3D viewer,
- Download generated protein structures.

<p align="center">
  <img src="images/BioEmu_Protein3D_Animated.gif" style="width: 70%; max-width: 800px; min-width: 300px;">
</p>

> [!NOTE]
> Provided companion Docker image is optimised for **CPU** use and doesn't require **GPU** compute.

## 📑 Table of contents:
- [Part 1: Local Use of Companion Docker Image](#part-1-local-use-of-companion-docker-image)
- [Part 2: Cloud Deployment to Azure Web App](#part-2-cloud-deployment-to-azure-web-app)
- [Part 3: User Experience - Gradio UI](#part-3-user-experience---gradio-ui)
- [Demo videos on YouTube](#demo-videos-on-youtube)
- [Acknowledgments and Citation](#acknowledgments-and-citation)

## Part 1: Local Use of Companion Docker Image
1. To pull the image locally, use the Docker command below. Alternatively, if you are deploying to an Azure Web app, skip this and move to Part 2.
``` PowerShell
docker pull ghcr.io/lazauk/bioemu-webapp:latest
```
2. You can start the Docker container with the following command:
``` PowerSehll
docker run -p 7860:7860 ghcr.io/lazauk/bioemu-webapp:latest
```
3. The Web app can now be accessed at http://127.0.0.1:7860.

> [!NOTE]
> Your machine should support the Docker Engine. On Windows, consider installing Docker Desktop.

## Part 2: Cloud Deployment to Azure Web App
1. Create a new Azure Web App and set the source container to:
    - `Image Source`: Other container registries
    - `Access type`: Public
    - `Registry server URL`: https://ghcr.io
    - `Image and tag`: lazauk/bioemu-webapp:latest
    - `Port`: 7860
![Web App - Container Config](images/BioEmu_WebApp_Config.png)
2. Because of the size of the Docker container, it may take a few minutes to pull the image during initial setup. You can verify the deployment status in *Deployment -> Deployment Center -> Logs* settings of your Azure Web app as shown below:
![Web App - Container Deploy](images/BioEmu_WebApp_Deploy.png)

> [!WARNING]
> For optimal performance, select an App Service Plan with at least 32GB of RAM.

## Part 3: User Experience - Gradio UI
1. Enter a protein sequence or choose from sample sequences;
2. Set the number of protein structures to predict;
3. Click "Generate Structures";
4. Examine the first 3D protein structure in the viewer;
5. Download the resulting structures if required.

## Demo videos on YouTube
🎥 You can see the solution in action on the following ["Predicting 3D Protein Structures with BioEmu AI"](https://youtu.be/k2yeQkbmGsg) YouTube video.

🎙️ You may also listen to the audio podcast, ["3D Protein Structures with BioEmu: Behind the Code"](https://youtu.be/szketeLILdc), where we review the mechanics of the provided Web UI wrapper in a beginner-friendly way.

## Acknowledgments and Citation
- [Microsoft BioEmu](https://github.com/microsoft/bioemu) - inference code and AI model weights;
- [Gradio](https://gradio.app/) - the Web UI framework.

The original BioEmu paper by the Microsoft Research team:
``` bibtex
@article {BioEmu2024,
    author = {Lewis, Sarah and Hempel, Tim and Jim{\'e}nez-Luna, Jos{\'e} and Gastegger, Michael and Xie, Yu and Foong, Andrew Y. K. and Satorras, Victor Garc{\'\i}a and Abdin, Osama and Veeling, Bastiaan S. and Zaporozhets, Iryna and Chen, Yaoyi and Yang, Soojung and Schneuing, Arne and Nigam, Jigyasa and Barbero, Federico and Stimper, Vincent and Campbell, Andrew and Yim, Jason and Lienen, Marten and Shi, Yu and Zheng, Shuxin and Schulz, Hannes and Munir, Usman and Clementi, Cecilia and No{\'e}, Frank},
    title = {Scalable emulation of protein equilibrium ensembles with generative deep learning},
    year = {2024},
    doi = {10.1101/2024.12.05.626885},
    journal = {bioRxiv}
}
```
