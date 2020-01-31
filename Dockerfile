FROM kbase/kbase:sdkbase2.latest
MAINTAINER KBase Developer

COPY ./ /kb/module

RUN mkdir -p /kb/module/work
RUN chmod -R 777 /kb/module
WORKDIR /kb/module

RUN cpanm --installdeps .

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

# CMD [ ]
