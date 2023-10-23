@ECHO OFF
TITLE BpresenT Launcher
REM Install and Launch script for Windows

REM Check if docker is available
docker version >nul
IF %ERRORLEVEL% NEQ 0 GOTO no_docker

REM Check if image is available
docker image inspect bpresent:latest >nul
IF %ERRORLEVEL% NEQ 0 GOTO no_image

ECHO ^>Starting BpresenT container
REM Start BpresentContainer, if it doesnt exist, pull image if needed, then run a new container
docker start BpresentContainer >nul || docker run --name BpresentContainer --pull missing -p 127.0.0.1::8000 -d bpresent:latest >nul
REM Give image/container 3 seconds to spin up
TIMEOUT /t 3 /nobreak >nul

ECHO ^>Opening in browser
REM Find bound port on host
@SET HOSTPORT=
FOR /F %%I IN ('docker inspect --format="{{(index (index .NetworkSettings.Ports \"8000/tcp\") 0).HostPort}}" BpresentContainer') DO @SET "HOSTPORT=%%I"
REM Open found address in default browser
RUNDLL32 url.dll,FileProtocolHandler http://127.0.0.1:%HOSTPORT%

TIMEOUT /t 1 /nobreak >nul
ECHO ^>[1;91mPress ANY key to stop running BpresenT container...[0m
REM Pause for user input, then stop container
PAUSE >nul
docker stop BpresentContainer >nul
GOTO:eof

:no_docker
ECHO ^>[1;91mDocker not found! Make sure it's running and in %%PATH%%.[0m
PAUSE
GOTO:eof

:no_image
ECHO ^>[1;91mBprenT Image not found![0m
PAUSE
GOTO:eof
