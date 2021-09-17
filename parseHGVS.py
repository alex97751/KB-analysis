import sys

import hgvs.parser
import hgvs.validator
import hgvs.dataproviders.uta
import hgvs.sequencevariant
import hgvs.edit
import hgvs.posedit


def setHGVS(variant):
    try:
        v = parseVariant(variant.hgvs)
        variant.ac = v.ac
        variant.posedit = v.type + "." + str(v.posedit)
        variant.hgvs_error = None
        variant.valid = validate(v, variant)
        # t = parseVariant(variant)
        # print(t.ac)
    except KeyboardInterrupt:
        print("Exiting.. \n")
        sys.exit(0)
    except Exception:
        # hgvs probably not valid
        variant.valid = False
        variant.hgvs_error = "NOT VALID HGVS"

    return


# return new hgvs variant
def parseVariant(variant):
    hp = hgvs.parser.Parser()
    var_g = hp.parse_hgvs_variant(variant)

    return var_g


# check if mutation on transcript is correct
def validate(var_g, variant):
    hdp = hgvs.dataproviders.uta.connect()
    vr = hgvs.validator.Validator(hdp)
    try:
        variant.hgvs_error = None
        return vr.validate(var_g)
        #print("is valid: " + str(vBool))
    except hgvs.exceptions.HGVSError as e:
        #print(e)
        variant.hgvs_error = e
        return False


# build a new hgvs variant from scratch
def newVariant(ac, type, ref, alt, base):
    edit = hgvs.edit.NARefAlt(ref='G', alt='A')
    start = hgvs.location.BaseOffsetPosition(base=1413, offset=-0, datum=hgvs.location.Datum.CDS_START)
    # end = hgvs.location.BaseOffsetPosition(base=22,datum=hgvs.location.Datum.CDS_END)
    iv = hgvs.location.Interval(start=start, end=start)
    var_g = hgvs.sequencevariant.SequenceVariant(ac='NM_014855.3', type='c', posedit=hgvs.posedit.PosEdit(pos=iv, edit=edit))

    print("New variant: " + str(var_g))
    return var_g

# DEBUG
# setHGVS("ENST00000369535:c.182delAAinsTG")  # same as ENST00000369535:c.182AA>TG
