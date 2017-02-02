FROM kingb12/fba_tools:latest

RUN mkdir /kb/module/lib/model-morphing
COPY ./ /kb/module/lib/model-morphing
RUN ["pip", "install", "-e", "/kb/module/lib/model-morphing"]

RUN ["chmod", "+x", "./lib/model-morphing/scripts/entrypoint.sh"]
RUN ["chmod", "+x", "./lib/model-morphing/scripts/run_async.sh"]
RUN ["chmod", "+x", "./lib/model-morphing/scripts/login.sh"]
RUN ["chmod", "+x", "./lib/model-morphing/all_the_models.sh"]

ENTRYPOINT [ "./lib/model-morphing/scripts/entrypoint.sh" ]
CMD ["login"]
