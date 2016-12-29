FROM kingb12/fba_tools:latest

RUN mkdir /kb/module/lib/model-morphing
COPY ./ /kb/module/lib/model-morphing

RUN ["pip", "install", "ipython"]
RUN pip install matplotlib
RUN pip install matplotlib_venn
RUN apt-get -y install vim

RUN ["chmod", "+x", "./lib/model-morphing/scripts/entrypoint.sh"]
RUN ["chmod", "+x", "./lib/model-morphing/scripts/run_async.sh"]
RUN ["chmod", "+x", "./lib/model-morphing/scripts/login.sh"]

ENTRYPOINT [ "./lib/model-morphing/scripts/entrypoint.sh" ]
