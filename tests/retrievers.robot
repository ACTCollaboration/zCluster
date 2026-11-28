*** Settings ***
Documentation       Retriever status checks: query each supported database for galaxies in a small box at
...                 a position inside its footprint, to find out which remote services are still working.
...
...                 These checks are meant to be run manually, on occasion, to find out if a service has been
...                 retired, moved, or changed - they are not part of the routine (regression) tests, and
...                 can take a long time to run. See docs/testing.rst.
...
...                 Nothing is cached between runs (each query uses a fresh, empty cache), so every test
...                 contacts the remote service. Tests are skipped if the credentials or Python modules a
...                 retriever needs are missing.
...
...                 Run from the repository root with:
...                 robot --outputdir tests/results tests/retrievers.robot
...
...                 Useful options:
...                 --include SDSSDR12              check a single database (tags are database names)
...                 --exclude slow                  skip the slow DECaLS tractor catalog downloads
...                 --exclude credentials           skip databases that need a login or token
Library             ZClusterLibrary.py    ${ZCLUSTER_ROOT}
Resource            common.resource
Test Template       Retriever Should Work
Test Timeout        15 minutes    # safety net - queries time out much sooner


*** Variables ***
${ZCLUSTER_ROOT}        ${CURDIR}/..
${HALF_BOX_DEG}         0.0667    # 4 arcmin
${MIN_GALAXIES}         10
${QUERY_TIMEOUT}        90        # seconds - a working service should respond well within this


*** Test Cases ***    DATABASE    RA_DEG    DEC_DEG    TIMEOUT
# Positions are 400SDSmall.csv clusters inside each survey footprint (cluster name in the comment),
# except for DP1, which contains no 400SD clusters. Rubin DP0 (simulated data) is not tested.
SDSS DR12             SDSSDR12    29.8258    0.5025    # 20
    [Tags]    SDSSDR12    SDSS

SDSS DR16             SDSSDR16    29.8258    0.5025    # 20
    [Tags]    SDSSDR16    SDSS

Pan-STARRS1 DR2       PS1         29.8258    0.5025    # 20
    [Tags]    PS1    credentials
    [Setup]    Skip Without CasJobs Login

DECaLS DR8            DECaLSDR8    29.8258    0.5025    600    # 20 (downloads whole tractor bricks)
    [Tags]    DECaLSDR8    DECaLS    slow

DECaLS DR9            DECaLSDR9    29.8258    0.5025    600    # 20 (downloads whole tractor bricks)
    [Tags]    DECaLSDR9    DECaLS    slow

DECaLS DR10           DECaLSDR10    29.8258    0.5025    600    # 20 (downloads whole tractor bricks)
    [Tags]    DECaLSDR10    DECaLS    slow

DECaLS DR10 (Data Lab)    DL_DECaLSDR10    29.8258    0.5025    # 20
    [Tags]    DL_DECaLSDR10    DECaLS    datalab

DECaLS DR11 (Data Lab)    DL_DECaLSDR11    29.8258    0.5025    # 20
    [Tags]    DL_DECaLSDR11    DECaLS    datalab

DELVE DR2 (Data Lab)    DELVEDR2    29.8258    0.5025    # 20
    [Tags]    DELVEDR2    datalab

DES DR1               DESDR1      25.3846    -30.5783    # 16
    [Tags]    DESDR1    DES    credentials
    [Setup]    Skip Without DES Login

DES DR2               DESDR2      25.3846    -30.5783    # 16
    [Tags]    DESDR2    DES    credentials
    [Setup]    Skip Without DES Login

KiDS DR4              KiDSDR4     20.6500    -28.5347    # 11
    [Tags]    KiDSDR4

ATLAS DR4             ATLASDR4    14.2292    -22.2308    # 6
    [Tags]    ATLASDR4

Rubin DP1             RubinDP1    53.1300    -28.1000    # ECDFS field centre (no 400SD clusters in DP1)
    [Tags]    RubinDP1    rubin    credentials
    [Setup]    Skip Without RSP Token

Rubin DP2             RubinDP2    7.4608    13.5039    # 1
    [Tags]    RubinDP2    rubin    credentials
    [Setup]    Skip Without RSP Token


*** Keywords ***
Retriever Should Work
    [Arguments]    ${database}    ${ra_deg}    ${dec_deg}    ${timeout}=${QUERY_TIMEOUT}
    ${n}=    Retriever Should Return Galaxies    ${database}    ${ra_deg}    ${dec_deg}    ${HALF_BOX_DEG}
    ...    ${OUTPUT DIR}/retrievers/${database}    ${MIN_GALAXIES}    ${timeout}
    Log    ${database}: ${n} galaxies    console=True
