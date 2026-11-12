"""Verify current saved results; no fitting, model solve or upstream reconstruction."""
from pathlib import Path
import csv,gzip,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'current_results'
def read(rel):
    p=DATA/rel
    opener=gzip.open if p.suffix=='.gz'else open
    with opener(p,'rt',encoding='utf-8-sig',newline='')as f:return list(csv.DictReader(f))
def main():
    inventory=json.loads((DATA/'manifest.json').read_text(encoding='utf-8'))
    for row in inventory['files']:
        p=ROOT/row['path']
        assert p.is_file()and hashlib.sha256(p.read_bytes()).hexdigest()==row['sha256'],row['path']
    response=read('data/response/environment_response.csv.gz')
    assert len(response)==67452,'Current response row identity differs'
    predictions=read('data/evaluation/formal_predictions_164.csv')
    assert len(predictions)==164 and len({r['sample_id']for r in predictions})==164
    reference=read('data/evaluation/s9_current_main.csv')
    output=[]
    for ref in reference:
        model=ref['model']
        errors=[math.log10(float(r[model])+.01)-math.log10(float(r['surface_LODhalf'])+.01)for r in predictions]
        mae=sum(abs(v)for v in errors)/len(errors);bias=sum(errors)/len(errors)
        maxpred=max(float(r[model])for r in predictions)
        assert len(errors)==int(ref['n'])==164
        assert abs(mae-float(ref['MAE']))<1e-12,(model,mae,ref['MAE'])
        assert abs(bias-float(ref['bias']))<1e-12,(model,bias,ref['bias'])
        assert abs(maxpred-float(ref['maximum_prediction_ugL']))<1e-10
        output.append({'model':model,'n':164,'MAE_log10_offset001':mae,'bias_log10_offset001':bias})
    print(json.dumps({'status':'CURRENT_SAVED_RESULT_RESCORING_PASS','response_rows':len(response),'metrics':output,'fits':0,'model_solves':0,'scope':'Archived current station-date predictions only; does not regenerate ecological states or upstream metabolic responses.'},indent=2))
if __name__=='__main__':main()
