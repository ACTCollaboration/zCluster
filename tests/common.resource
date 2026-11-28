*** Settings ***
Documentation       Keywords shared by the zCluster test suites.


*** Keywords ***
Skip Without RSP Token
    Skip If    not os.path.exists(os.path.expanduser('~/.rsp-tap.token'))
    ...    No RSP token in ~/.rsp-tap.token - see https://rsp.lsst.io/guides/auth/creating-user-tokens.html

Skip Without CasJobs Login
    Skip If    not (os.environ.get('CASJOBS_WSID') and os.environ.get('CASJOBS_PW'))
    ...    CASJOBS_WSID and CASJOBS_PW not set - see http://ps1images.stsci.edu/ps1_dr2_query.html

Skip Without DES Login
    Skip If    not os.path.exists(os.path.expanduser('~/.desservices.ini'))
    ...    No DES database login in ~/.desservices.ini - see https://github.com/des-labs/easyaccess
