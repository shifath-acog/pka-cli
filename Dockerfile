FROM dockerhub.aganitha.ai:4443/chem/pka_pipeline:ovh4

CMD ["jupyter" "lab" "--ip=0.0.0.0" "--no-browser" "--port=8888" "--NotebookApp.token=''" "--NotebookApp.password=''" "--allow-root" "--ContentsManager.allow_hidden=True"]
