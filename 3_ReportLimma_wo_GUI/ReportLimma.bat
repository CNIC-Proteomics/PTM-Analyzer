@echo off
setlocal

:: Get the path of the current batch script
set "base_path=%~dp0"

:: Initialize variables for each expected argument
set "input="
set "samples="
set "config="
set "output="
set "show_help=false"


goto parse_args



:: Function to print help
:print_help
echo Usage: %~nx0 -i ^<input_file^> -s ^<samples_file^> -c ^<config_file^> -o ^<output_file^>
echo.
echo Options:
echo     -i    Path to the isanxot report file
echo     -s    Path to the sample table file
echo     -c    Path to the YAML configuration file
echo     -o    Path to the output file
echo     -h, --help    Show this help message and exit
echo.
goto end_script 1



:: Parse arguments
:parse_args
if "%~1"=="" goto end_parse
if "%~1"=="-i" (
    set "input=%~2"
    shift
) else if "%~1"=="-s" (
    set "samples=%~2"
    shift
) else if "%~1"=="-c" (
    set "config=%~2"
    shift
) else if "%~1"=="-o" (
    set "output=%~2"
    shift
) else if "%~1"=="-h" (
    set "show_help=true"
) else if "%~1"=="--help" (
    set "show_help=true"
)
shift
goto parse_args

:end_parse



:: Show help if requested
if "%show_help%"=="true" (
    goto print_help
)


:: Validate required parameters
if "%input%"=="" (
    echo Missing '-i' input argument
    goto print_help
)
if "%samples%"=="" (
    echo Missing '-s' samples argument
    goto print_help
)
if "%config%"=="" (
    echo Missing '-c' config argument
    goto print_help
)
if "%output%"=="" (
    echo Missing -o output argument
    goto print_help
)



:: Set paths
set "Rexe_relative=R-Portable\App\R-Portable\bin\x64\Rscript.exe"
set "Rexe_path=%base_path%%Rexe_relative%"
set "Rscript=%base_path%app_wo_GUI.R"



:: Execute command
echo.
echo ** Executing command...
echo "%Rexe_path%" --vanilla "%Rscript%" -i "%input%" -s "%samples%" -c "%config%" -o "%output%"
"%Rexe_path%" --vanilla "%Rscript%" -i "%input%" -s "%samples%" -c "%config%" -o "%output%"

goto end_script



:end_script
endlocal
exit /b 0
