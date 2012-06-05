from sqlalchemy import DateTime, MetaData, Column, Table, ForeignKey, Integer, String, Sequence, Boolean, func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import datetime
import json

Base = declarative_base()

class Releases_Table(Base):
    __tablename__ = 'releases'

    id = Column(Integer, Sequence('releases_id_seq'), primary_key=True, nullable=False)
    category = Column(String, nullable=False)
    subcategory = Column(String, nullable=False)
    release_name = Column(String, nullable=False)
    version = Column(Integer, nullable=False)
    date = Column(DateTime, nullable=False)
    status_kind = Column(String, nullable=False)

    def __init__(self, category, subcategory, release_name, version, date, status_kind):
        self.category = category
        self.subcategory = subcategory
        self.release_name = release_name
        self.version = version
        self.date = date
        self.status_kind = status_kind

    def __repr__(self):
       return "<Releases_Table('%s','%s', '%s', '%s', '%s', '%s')>" % (self.category, self.subcategory, self.release_name, self.version, self.date, self.status_kind)
       
class Releases_LV_Table(Base):
    __tablename__ = 'releases_lv'

    id = Column(Integer, primary_key=True, nullable=False)
    category = Column(String, nullable=False)
    subcategory = Column(String, nullable=False)
    release_name = Column(String, nullable=False)
    version = Column(Integer, nullable=False)
    date = Column(DateTime, nullable=False)
    status_kind = Column(String, nullable=False)

    def __init__(self, id, category, subcategory, release_name, version, date, status_kind):
        self.id = id
        self.category = category
        self.subcategory = subcategory
        self.release_name = release_name
        self.version = version
        self.date = date
        self.status_kind = status_kind

    def __repr__(self):
       return "<Releases_LV_Table('%s','%s', '%s', '%s', '%s', '%s')>" % (self.category, self.subcategory, self.release_name, self.version, self.date, self.status_kind)
       
class Status_Table(Base):
    __tablename__ = 'status'
    
    id = Column(Integer, ForeignKey("releases.id"), primary_key=True, nullable=False)
    validation_status = Column(String, default="NOT YET DONE", nullable=False)
    comments = Column(String)
    links = Column(String)
    user_name = Column(String, nullable=False)

    def __init__(self, id, validation_status, comments, links, user_name):
        self.id = id
        self.validation_status = validation_status
        self.comments = comments
        self.links = links
        self.user_name = user_name

    def __repr__(self):
       return "<Status_Table('%d','%s', '%s', '%s', '%s')>" % (self.id, self.validation_status, self.comments, self.links, self.user_name)
       
class Status_LV_Table(Base):
    __tablename__ = 'status_lv'
    
    id = Column(Integer, ForeignKey("releases_lv.id"), primary_key=True, nullable=False)
    validation_status = Column(String, default="NOT YET DONE", nullable=False)
    comments = Column(String)
    links = Column(String)
    user_name = Column(String, nullable=False)

    def __init__(self, id, validation_status, comments, links, user_name):
        self.id = id
        self.validation_status = validation_status
        self.comments = comments
        self.links = links
        self.user_name = user_name

    def __repr__(self):
       return "<Status_LV_Table('%d','%s', '%s', '%s', '%s')>" % (self.id, self.validation_status, self.comments, self.links, self.user_name)

class Users_Table(Base):
    __tablename__ = 'users'
    
    user_name = Column(String, nullable=False, primary_key=True)
    email = Column(String)
    admin = Column(Boolean, default=False)
    validator = Column(Boolean, default=False)

    def __init__(self, user_name, email, admin, validator):
        self.user_name = user_name
        self.email = email
        self.admin = admin
        self.validator = validator

    def __repr__(self):
       return "<Users_Table('%s', '%s', '%b', '%b')>" % (self.user_name, self.email, self.admin, self.validator)
       
class User_Rights_Table(Base):
    __tablename__ = 'user_rights'
    
    id = Column(Integer, Sequence('id_seq'), primary_key=True, nullable=False)
    user_name = Column(String, ForeignKey("users.user_name"), nullable=False)
    category = Column(String, nullable=False)
    subcategory = Column(String, nullable=False)
    status_kind = Column(String, nullable=False)

    def __init__(self, user_name, category, subcategory, status_kind):
        self.user_name = user_name
        self.category = category
        self.subcategory = subcategory
        self.status_kind = status_kind

    def __repr__(self):
       return "<User_Rights_Table('%s','%s', '%s', '%s')>" % (self.user_name, self.category, self.subcategory, self.status_kind)

RELEASE_NAME = "RELEASE_NAME"
VALIDATION_STATUS = "VALIDATION_STATUS"
COMMENTS = "COMMENTS"
LINKS = "LINKS"
META_DATE = "META_DATE"
USER_NAME = "USER_NAME"
CATEGORY = "CATEBORY"
SUBCATEGORY = "SUBCATEGORY"
    
possible_status_list = ["OK", 
                        "NOT YET DONE", 
                        "FAILURE", 
                        "CHANGES EXPECTED", 
                        "IN PROGRESS", 
                        "OK TO BE SIGNED-OFF BY THE VALIDATORS", 
                        "FAILURE TO BE SIGNED-OFF BY THE VALIDATORS", 
                        None]
possible_category_list = ["Reconstruction", "HLT", "PAGs"]
possible_subcatrgory_list = ["Data", "FullSim", "FastSim"] 

reconstruction_status_list = ["CSC", "TAU", "TRACKING", "BTAG", "JET", "ECAL", "RPC", "PHOTON", "MUON", "MET", "ELECTRON", "TK", "HCAL", "DT", "SUMMARY"]
hlt_status_list = ["TAU", "JET", "HIGGS", "TOP", "MUON", "PHOTON", "MET", "ELECTRON", "EXOTICA", "SUSY", "HEAVYFLAVOR", "SUMMARY"]
pags_status_list = ["B", "HIGGS", "FWD", "TOP", "SMP", "EXOTICA", "SUSY", "SUMMARY"]

# Returns validation statuses in JSON key-value form, found by release category, subcategory and name 
def getReleaseShortInfo(cat, sub_cat, rel_name, Session):
    session = Session()
    try:
        info_dict = {}
        for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                filter(Releases_Table.subcategory == sub_cat).\
                                                filter(Releases_Table.release_name == rel_name):
            for j in session.query(Status_Table).filter(Status_Table.id == i.id):
                info_dict[i.status_kind] = j.validation_status
        session.close()
        return json.dumps(info_dict)
    except Exception as e:
        session.close()
        print e

# Returns validation status by column of release
def getStatus(cat, sub_cat, rel_name, status_kind, Session):
    session = Session()
    try:
        id = None
        for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                filter(Releases_Table.subcategory == sub_cat).\
                                                filter(Releases_Table.release_name == rel_name).\
                                                filter(Releases_Table.status_kind == status_kind):
            id = i.id
        for i in session.query(Status_Table).filter(Status_Table.id == id):
            session.close()
            return i.validation_status
        session.close()
        return "Error! Unknown status kind"
    except Exception as e:
        session.close()
        print e
    
# Returns release status details in JSON key-value form, found by release category, subcategory and name and kind of validation status
def getReleaseDetails(cat, sub_cat, rel_name, status_kind, Session):
    session = Session()
    try:
        version = 0
        info_dict = {}
        for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                filter(Releases_Table.subcategory == sub_cat).\
                                                filter(Releases_Table.release_name == rel_name).\
                                                filter(Releases_Table.status_kind == status_kind):
            for j in session.query(Status_Table).filter(Status_Table.id == i.id):
                info_dict[RELEASE_NAME] = i.release_name
                info_dict[VALIDATION_STATUS] = j.validation_status
                info_dict[COMMENTS] = j.comments
                info_dict[LINKS] = j.links
                info_dict[META_DATE] = i.date.strftime("%Y-%m-%d %H:%M:%S")
                info_dict[USER_NAME] = j.user_name
        session.close()
        return json.dumps(info_dict)
    except Exception as e:
        session.close()
        print e

# Returns release name(s) in JSON list, found by keyword
def search(regexp, Session):
    session = Session()
    try:
        regexp = regexp.lower()
        release_list = []
        new_list = []
        for i in session.query(Releases_Table.release_name, func.count(Releases_Table.release_name)).group_by(Releases_Table.release_name).all():
            if i.release_name not in release_list:
                release_list.append(i.release_name)
        if regexp == "*" or regexp == "" or regexp == None:
            release_list.sort(key=lambda x: x.lower())
            session.close()
            return json.dumps(release_list)
        elif not "*" in regexp:
            for name in release_list:
                if regexp.lower() in name.lower():
                    new_list.append(name)
            new_list.sort(key=lambda x: x.lower())
            session.close()
            return json.dumps(new_list)
        else:
            session.close()
            return json.dumps([])
    except Exception as e:
        session.close()
        print e
 
# Returns all history of releases
def getAllHistory(Session):
    session = Session()
    try:
        release_list = []
        for i in session.query(Releases_Table.release_name, func.count(Releases_Table.release_name)).group_by(Releases_Table.release_name).all():
                if i.release_name not in release_list:
                    release_list.append(i.release_name)
        session.close()    
        return getHistory(release_list, possible_subcatrgory_list, reconstruction_status_list, hlt_status_list, pags_status_list, Session)
    except Exception as e:
        session.close()
        print e
        return "Error in database reading!"
   
# Returns history of releases
def getHistory(release_list, subcategory_list, reco_list, hlt_list, pags_list, Session):
    session = Session()
    try:
        history_dict = {}
        for rel_name in release_list:
            cat_dict = {}
            for i in possible_category_list:
                sub_cat_dict = {}
                for j in possible_subcatrgory_list:
                    sub_cat_dict[j] = rel_name
                cat_dict[i] = sub_cat_dict
            A = {}
            for i in cat_dict:
                B = {}
                for j in cat_dict[i]:
                    C = {}
                    if j not in subcategory_list:
                        continue
                    for k in session.query(Releases_LV_Table.version, func.count(Releases_LV_Table.version)).filter(Releases_LV_Table.release_name == rel_name).\
                                                            filter(Releases_LV_Table.category == i).\
                                                            filter(Releases_LV_Table.subcategory == j).\
                                                            group_by(Releases_LV_Table.version).all():
                        D = {}
                        for l in session.query(Releases_LV_Table).filter(Releases_LV_Table.release_name == rel_name).\
                                                                filter(Releases_LV_Table.category == i).\
                                                                filter(Releases_LV_Table.subcategory == j).\
                                                                filter(Releases_LV_Table.version == k.version):
                            E = {}
                            if i == possible_category_list[0]: # Reconstruction
                                if l.status_kind not in reco_list:
                                    continue
                            if i == possible_category_list[1]: # HLT
                                if l.status_kind not in hlt_list:
                                    continue
                            if i == possible_category_list[2]: # PAGs
                                if l.status_kind not in pags_list:
                                    continue
                            for m in session.query(Status_LV_Table).filter(Status_LV_Table.id == l.id):
                                E[VALIDATION_STATUS] = m.validation_status
                                E[META_DATE] = l.date.strftime("%Y-%m-%d %H:%M:%S")
                                E[COMMENTS] = m.comments
                                E[LINKS] = m.links
                                E[USER_NAME] = m.user_name
                            if E: D[l.status_kind] = E
                        if D: C[k.version] = D
                    if C: B[j] = C
                if B: A[i] = B
            history_dict[rel_name] = A
        session.close()
        return json.dumps(history_dict)
    except Exception as e:
        session.close()
        print e
       
# Creates new release. Parameters:   cat - release category
#                                    sub_cat - release subcategory
#                                    rel_name - release neme
#                                    list - JSON list contains status information. For example: {SMP:{VALIDATION STATUS:"OK", COMMENTS:"...", LINKS:"...", USER_NAME:"..."}}
def newRelease(cat, sub_cat, rel_name, dict_json, Session, *args):
    session = Session()
    try:
        version = None
        if len(args) > 0:
            version = args[0]
            session = args[1]
        else:
            version = 1
            session = Session()
        if cat not in possible_category_list or sub_cat not in possible_subcatrgory_list:
            session.close()
            return "Error! Wrong category or subcategory"
        
        for i in session.query(Releases_Table).filter(Releases_Table.category == cat).filter(Releases_Table.subcategory == sub_cat).filter(Releases_Table.release_name == rel_name):
            if i != None:
                session.close()
                return "Error! Release already exists"
        dict = json.loads(dict_json)
        list = dict.keys()
        date = datetime.datetime.now()
        
        for status_kind in list:
            if dict[status_kind][VALIDATION_STATUS] not in possible_status_list:
                session.close()
                return "Error! Illegal validation status"
            release = Releases_Table(cat, sub_cat, rel_name, version, date, status_kind)
            session.add(release)
            id = None
            for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                    filter(Releases_Table.subcategory == sub_cat).\
                                                    filter(Releases_Table.release_name == rel_name).\
                                                    filter(Releases_Table.status_kind == status_kind):
                id = i.id
            user_name = dict[status_kind][USER_NAME]
            if user_name == "":
                user_name = "Unknown name"
            status = Status_Table(id, dict[status_kind][VALIDATION_STATUS], dict[status_kind][COMMENTS], dict[status_kind][LINKS], user_name)
            session.add(status)
            release_lv = Releases_LV_Table(id, cat, sub_cat, rel_name, version, date, status_kind)
            session.add(release_lv)
            id_lv = None
            for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                    filter(Releases_Table.subcategory == sub_cat).\
                                                    filter(Releases_Table.release_name == rel_name).\
                                                    filter(Releases_LV_Table.status_kind == status_kind):
                id_lv = i.id
            status_lv = Status_LV_Table(id, dict[status_kind][VALIDATION_STATUS], dict[status_kind][COMMENTS], dict[status_kind][LINKS], user_name)
            session.add(status_lv)
        if len(args) > 0:
            session.commit()
            session.close()
            return "True", session
        else:
            session.commit()
            session.close()
            return "True"
    except Exception as e:
        session.close()
        print e

# Changes validation status of given release
def changeStatus(cat, sub_cat, rel_name, status_kind, new_status, new_comment, new_user_name, new_links, Session):
    session = Session()
    try:
        date = datetime.datetime.now()
        version = None 
        for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                filter(Releases_Table.subcategory == sub_cat).\
                                                filter(Releases_Table.release_name == rel_name):
            version = i.version
        if version == None:
            session.close()
            return "Error! Release not exists"
        version = version + 1
        dict = {}
        status_id_for_delete = []
        for i in session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                                filter(Releases_Table.subcategory == sub_cat).\
                                                filter(Releases_Table.release_name == rel_name):
            status_dict = {}
            status_id_for_delete.append(i.id)
            for j in session.query(Status_Table).filter(Status_Table.id == i.id):
                status_dict[VALIDATION_STATUS] = j.validation_status
                status_dict[COMMENTS] = j.comments
                status_dict[USER_NAME] = j.user_name
                status_dict[LINKS] = j.links
            dict[i.status_kind] = status_dict
        dict[status_kind][VALIDATION_STATUS] = new_status
        dict[status_kind][COMMENTS] = new_comment
        dict[status_kind][LINKS] = new_links
        dict[status_kind][USER_NAME] = new_user_name
        for i in status_id_for_delete:
            session.query(Status_Table).filter(Status_Table.id == i).delete()
        session.query(Releases_Table).filter(Releases_Table.category == cat).\
                                        filter(Releases_Table.subcategory == sub_cat).\
                                        filter(Releases_Table.release_name == rel_name).delete()        
        response = newRelease(cat, sub_cat, rel_name, json.dumps(dict), Session, version, session)
        if response[0] == "True":
            response[1].commit()
            session.close()
            return "True"
        else:
            session.rollback()
            session.close()
            return "Error!"
    except Exception as e:
        session.close()
        print e
  
#=======================USERS===========================

# Adds user to database
def addUser(user_name, email, Session):
    session = Session()
    try:
        user_name = user_name.lower()
        exists = False
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            if i.user_name == user_name:
                exists = True
        if exists:
            session.close()
            return "Error! Already exists"
        user = Users_Table(user_name, email, None, None)
        session.add(user)
        session.commit()
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            if i.user_name == user_name:
                exists = True
        if exists:
            session.close()
            return "True"
        else:
            session.close()
            return "False"
    except Exception as e:
        session.close()
        print e

# Removes user from database
def removeUser(user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        exists = False
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            if i.user_name == user_name:
                exists = True
        if not exists:
            session.close()
            return "User not exist"
        session.query(User_Rights_Table).filter(User_Rights_Table.user_name == user_name).delete()
        session.query(Users_Table).filter(Users_Table.user_name == user_name).delete()
        session.commit()
        exists = False
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            if i.user_name == user_name:
                exists = True
        if not exists:
            session.close()
            return "True"
        else:
            session.close()
            return "False"
    except Exception as e:
        session.close()
        print e

# Returns emlail by username
def getUserEmail(user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            session.close()
            return i.email
    except Exception as e:
        session.close()
        print e
        
# Grants administrator rights
def grantAdminRights(user_name, Session):
    session = Session()
    user_name = user_name.lower()
    try:
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            i.admin = True
        session.commit()
        session.close()
        return "True"
    except Exception as e:
        session.close()
        print e
        return "False"

# Grants validator rights
def grantValidatorRights(user_name, Session):
    session = Session()
    user_name = user_name.lower()
    try:
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            i.validator = True
        session.commit()
        session.close()
        return "True"
    except Exception as e:
        session.close()
        print e
        return "False"

# Withdraws administrator rights
def withdrawAdminRights(user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            i.admin = False
        session.commit()
        session.close()
        return "True"
    except Exception as e:
        session.close()
        print e
        return "False"

# Withdraws validator rights
def withdrawValidatorRights(user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        for i in session.query(Users_Table).filter(Users_Table.user_name == user_name):
            i.validator = False
        session.commit()
        session.close()
        return "True"
    except Exception as e:
        session.close()
        print e
        return "False"

# Returns TRUE if user is administrator and FALSE - vice versa
def checkAdmin(user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        for i in session.query(Users_Table.admin).filter(Users_Table.user_name == user_name):
            session.close()
            return i[0]
        session.close()
        return False
    except Exception as e:
        session.close()
        print e
        return False

# Returns TRUE if user is validator and FALSE - vice versa
def checkValidator(user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        for i in session.query(Users_Table.validator).filter(Users_Table.user_name == user_name):
            session.close()
            return i[0]
        session.close()
        return False
    except Exception as e:
        session.close()
        print e
        return False
    
# Checks validator rights to modify release status
def checkValidatorRights(cat, sub_cat, status_kind, user_name, Session):
    user_name = user_name.lower()
    session = Session()
    try:
        for i in session.query(User_Rights_Table).filter(User_Rights_Table.category == cat).\
                                                    filter(User_Rights_Table.subcategory == sub_cat).\
                                                    filter(User_Rights_Table.user_name == user_name).\
                                                    filter(User_Rights_Table.status_kind == status_kind):
            if i.user_name == user_name:
                session.close()
                return True
        session.close()
        return "False"
    except Exception as e:
        session.close()
        print e

    
# Grants validator rights to modify release status    
def grantValidatorRightsForStatusKind(cat, sub_cat, status_kind, user_name, Session):
    session = Session()
    user_name = user_name.lower()
    try:
        validator = User_Rights_Table(user_name, cat, sub_cat, status_kind)
        session.add(validator)
        session.commit()
        session.close()
        return "True"
    except Exception as e:
        session.close()
        print e
        return "False"


# Grants validator rights to modify validation data
def grantValidatorRightsForStatusKindList(user_name, rec_dat_list, rec_ful_list, rec_fas_list, hlt_dat_list, hlt_ful_list, hlt_fas_list, pag_dat_list, pag_ful_list, pag_fas_list, Session):
    user_name = user_name.lower()
    flag = True

    for status_kind in rec_dat_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[0], possible_subcatrgory_list[0], status_kind, user_name, Session) # Reconstruction-Data
    if not flag:
        return "False"

    for status_kind in rec_ful_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[0], possible_subcatrgory_list[1], status_kind, user_name, Session) # Reconstruction-FullSim
    if not flag:
        return "False"

    for status_kind in rec_fas_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[0], possible_subcatrgory_list[2], status_kind, user_name, Session) # Reconstruction-FastSim
    if not flag:
        return "False"
    #----------------------
    for status_kind in hlt_dat_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[1], possible_subcatrgory_list[0], status_kind, user_name, Session) # HLT-Data
    if not flag:
        return "False"

    for status_kind in hlt_ful_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[1], possible_subcatrgory_list[1], status_kind, user_name, Session) # HLT-FullSim
    if not flag:
        return "False"

    for status_kind in hlt_fas_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[1], possible_subcatrgory_list[2], status_kind, user_name, Session) # HLT-FastSim
    if not flag:
        return "False"
    #----------------------
    for status_kind in pag_dat_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[2], possible_subcatrgory_list[0], status_kind, user_name, Session) # PAGs-Data
    if not flag:
        return "False"

    for status_kind in pag_ful_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[2], possible_subcatrgory_list[1], status_kind, user_name, Session) # PAGs-FullSim
    if not flag:
        return "False"

    for status_kind in pag_fas_list:
        if status_kind == "":
            continue
        flag = grantValidatorRightsForStatusKind(possible_category_list[2], possible_subcatrgory_list[2], status_kind, user_name, Session) # PAGs-FastSim
    if not flag:
        return "False"
    
    if flag:
        return "True"
    else:
        return "False"

# Returns users name(s) in JSON list, found by keyword
def searchUsers(regexp, Session):
    session = Session()
    try:
        regexp = regexp.lower()
        user_list = []
        new_list = []
        for i in session.query(Users_Table):
            if i.user_name not in user_list:
                user_list.append(i.user_name)
        if regexp == "*" or regexp == "" or regexp == None:
            user_list.sort(key=lambda x: x.lower())
            session.close()
            return json.dumps(user_list)
        elif not "*" in regexp:
            for name in user_list:
                if regexp.lower() in name.lower():
                    new_list.append(name)
            new_list.sort(key=lambda x: x.lower())
            session.close()
            return json.dumps(new_list)
        else:
            session.close()
            return json.dumps([])
    except Exception as e:
        session.close()
        print e

# Returns all information (in JSON) about users rights
def getAllUsersInfo(regexp, Session):
    user_list = json.loads(searchUsers(regexp, Session))
    session = Session()
    try:
        dict = {}
        admin_list = []
        validator_dict = {}
        for i in session.query(Users_Table):
            if i.admin == True and i.user_name in user_list:
                admin_list.append(i.user_name)
            elif i.validator == True and i.user_name in user_list:
                validator_dict[i.user_name] = None
        dict["admins"] = admin_list
        dict["validators"] = validator_dict
        
        validators = dict["validators"].keys()
        for validator in validators:
            dict["validators"][validator] = {"Reconstruction":{"Data":[], "FullSim":[], "FastSim":[]},
                                            "HLT":{"Data":[], "FullSim":[], "FastSim":[]},
                                            "PAGs":{"Data":[], "FullSim":[], "FastSim":[]}}
        for validator in validators:
            for row in session.query(User_Rights_Table).filter(User_Rights_Table.user_name == validator):
                dict["validators"][validator][row.category][row.subcategory].append(row.status_kind)

        session.close()
        return json.dumps(dict)
    except Exception as e:
        session.close()
        print e
        return "Error in database reading!"
