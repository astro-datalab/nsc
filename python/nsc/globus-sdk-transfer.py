#!/usr/bin/env python 

# Imports
import argparse
import globus_sdk
from globus_sdk.scopes import TransferScopes
from globus_sdk.tokenstorage import SimpleJSONFileAdapter
from globus_sdk import (
    AuthClient,
    TransferClient,
    ConfidentialAppAuthClient,
    RefreshTokenAuthorizer,)
import numpy as np
import os

if __name__=="__main__":

    my_file_adapter = SimpleJSONFileAdapter(os.path.expanduser("~/mytokens.json"))
    my_file_adapterC = SimpleJSONFileAdapter(os.path.expanduser("~/mytokensC.json"))
    
    parser = argparse.ArgumentParser()
    parser.add_argument("source_id") # name of file source location
    parser.add_argument("dest_id") # name of file destination location
    parser.add_argument("file_sourcebase") # basedir of file sources
    parser.add_argument("file_destbase") # basedir of file dests
    parser.add_argument("file_sources") # files to transfer
    parser.add_argument("file_dests")
    args = parser.parse_args()
    print("file_sources = ",args.file_sources.split(","))
    print("file_destss = ",args.file_dests.split(","))

    for fls,fld in zip(args.file_sources.split(','),args.file_dests.split(',')):
        print(args.file_sourcebase+fls,args.file_destbase+fld)

    hosts = np.array(["tempest","blackmore"])
    ids = np.array(["0dc1297f-9868-4c68-8637-c9b6bd65d3aa","5485832e-723e-4b52-8472-0410e90902ad"])

    source_id = ids[hosts==args.source_id][0]
    dest_id = ids[hosts==args.dest_id][0]

    #CLIENT_ID = "645d1f13-700b-41b9-930f-0b3e73f3357d" 
    CLIENT_ID = "ee825967-32c6-4906-9278-ace1bcbadfa2"
    auth_client = globus_sdk.NativeAppAuthClient(CLIENT_ID)

    # we will need to do the login flow potentially twice, so define it as a
    # function
    #
    # we default to using the Transfer "all" scope, but it is settable here
    # look at the ConsentRequired handler below for how this is used

    def login_and_get_transfer_client(*, scopes=TransferScopes.all):
        # note that 'requested_scopes' can be a single scope or a list
        # this did not matter in previous examples but will be leveraged in
        # this one
        
        if not my_file_adapter.file_exists():
            auth_client.oauth2_start_flow(requested_scopes=scopes,refresh_tokens=True)
            authorize_url = auth_client.oauth2_get_authorize_url()
            print(f"Please go to this URL and login:\n\n{authorize_url}\n")
    
            auth_code = input("Please enter the code here: ").strip()
            tokens = auth_client.oauth2_exchange_code_for_tokens(auth_code)
            transfer_tokens = tokens.by_resource_server["transfer.api.globus.org"]
            
            transfer_rt = transfer_tokens["refresh_token"]
            transfer_at = transfer_tokens["access_token"]
            expires_at_s = transfer_tokens["expires_at_seconds"]
            
            my_file_adapter.store(tokens)
            
        else:
            transfer_tokens = my_file_adapter.get_token_data("transfer.api.globus.org")
            transfer_rt = transfer_tokens["refresh_token"]
            transfer_at = transfer_tokens["access_token"]
            expires_at_s = transfer_tokens["expires_at_seconds"]
        
        
        
        transfer_authorizer = RefreshTokenAuthorizer(
        transfer_rt, auth_client, access_token=transfer_at, expires_at=expires_at_s
        )   
        
        # return the TransferClient object, as the result of doing a login
        return globus_sdk.TransferClient(
            authorizer=transfer_authorizer)
            
            
            
            
            
    def login_and_get_transfer_clientC(*, scopes=TransferScopes.all):
        # note that 'requested_scopes' can be a single scope or a list
        # this did not matter in previous examples but will be leveraged in
        # this one
        
        if not my_file_adapterC.file_exists():
            auth_client.oauth2_start_flow(requested_scopes=scopes,refresh_tokens=True)
            authorize_url = auth_client.oauth2_get_authorize_url()
            print(f"Please go to this URL and login:\n\n{authorize_url}\n")
    
            auth_code = input("Please enter the code here: ").strip()
            tokensC = auth_client.oauth2_exchange_code_for_tokens(auth_code)
            transfer_tokensC = tokensC.by_resource_server["transfer.api.globus.org"]
            
            transfer_rt = transfer_tokensC["refresh_token"]
            transfer_at = transfer_tokensC["access_token"]
            expires_at_s = transfer_tokensC["expires_at_seconds"]
            
            my_file_adapterC.store(tokensC)
            
        else:
            transfer_tokensC = my_file_adapterC.get_token_data("transfer.api.globus.org")
            transfer_rt = transfer_tokensC["refresh_token"]
            transfer_at = transfer_tokensC["access_token"]
            expires_at_s = transfer_tokensC["expires_at_seconds"]
        
        
        
        transfer_authorizer = RefreshTokenAuthorizer(
        transfer_rt, auth_client, access_token=transfer_at, expires_at=expires_at_s
        )   
        
        # return the TransferClient object, as the result of doing a login
        return globus_sdk.TransferClient(
            authorizer=transfer_authorizer)        
    
    
    
            
            
            
    
    
    # get an initial client to try with, which requires a login flow
    transfer_client = login_and_get_transfer_client()
    
    # now, try an ls on the source and destination to see if ConsentRequired
    # errors are raised
    consent_required_scopes = []
    
    
    def check_for_consent_required(target):
        try:
            transfer_client.operation_ls(target, path="/")
        # catch all errors and discard those other than ConsentRequired
        # e.g. ignore PermissionDenied errors as not relevant
        except globus_sdk.TransferAPIError as err:
            if err.info.consent_required:
                consent_required_scopes.extend(err.info.consent_required.required_scopes)
    
    
    check_for_consent_required(source_id)
    check_for_consent_required(dest_id)
    
    # the block above may or may not populate this list
    # but if it does, handle ConsentRequired with a new login
    if consent_required_scopes:
        print(
            "One of your endpoints requires consent in order to be used.\n"
            "You must login a second time to grant consents.\n\n"
        )
        transfer_client = login_and_get_transfer_clientC(scopes=consent_required_scopes)
    
    # from this point onwards, the example is exactly the same as the reactive
    # case, including the behavior to retry on ConsentRequiredErrors. This is
    # not obvious, but there are cases in which it is necessary -- for example,
    # if a user consents at the start, but the process of building task_data is
    # slow, they could revoke their consent before the submission step
    #
    # in the common case, a single submission with no retry would suffice
    
    task_data = globus_sdk.TransferData(
        source_endpoint=source_id, destination_endpoint=dest_id
    )
    for fls,fld in zip(args.file_sources.split(','),args.file_dests.split(',')):
        task_data.add_item(
             args.file_sourcebase+fls,
             args.file_destbase+fld,
        )
    #task_data.add_item(
    #    "/home/x25h971/nsc/instcal/v4/phot_tempest.py",
    #    "/phyx-nidever/kfas/nsc_meas/phot_tempest.py",
    #)
    #task_data.add_item(
    #    "/home/d81v711/helloWorld.class",  # source
    #    "/uit-rci/Nitasha-Test/helloWorld.class",  # dest
    #)
    
    
    def do_submit(client):
        task_doc = client.submit_transfer(task_data)
        task_id = task_doc["task_id"]
        print(f"submitted transfer, task_id={task_id}")
        return(task_id)
    
    
    try:
        task_id = do_submit(transfer_client)
    except globus_sdk.TransferAPIError as err:
        if not err.info.consent_required:
            raise
        print(
            "Encountered a ConsentRequired error.\n"
            "You must login a second time to grant consents.\n\n"
        )
        transfer_client = login_and_get_transfer_clientC(
            scopes=err.info.consent_required.required_scopes
        )
        task_id = do_submit(transfer_client)
        
    done = transfer_client.task_wait(task_id, timeout=600, polling_interval=60)
    if not done:
        print(f"{task_id} didn't successfully terminate!")
    else:
        print(f"{task_id} completed")
